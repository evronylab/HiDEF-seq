#!/usr/bin/env python3
"""Content identities and atomic, immutable publication of prepared artifacts.

This helper is deliberately independent of Nextflow's filename-based storeDir.
Consumers must use the returned directory, and must never discover products by
looking in a partially populated directory. See docs/workflow-optimization.md.
"""

import argparse
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


SCHEMA = 1
MANIFEST = "manifest.complete.json"


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("utf-8")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identify(spec):
    """Hash each relevant input once; never identify data solely by its name."""
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", spec["namespace"]):
        raise ValueError("namespace must be a simple, nonempty directory name")
    if spec["namespace"] in (".", ".."):
        raise ValueError("invalid namespace")
    unknown = set(spec) - {"namespace", "settings", "inputs", "scripts", "tools"}
    if unknown:
        raise ValueError("unknown identity fields: " + ", ".join(sorted(unknown)))
    identity = {"schema": SCHEMA, "namespace": spec["namespace"],
                "settings": json.loads(canonical(spec.get("settings", {}))),
                "tools": json.loads(canonical(spec.get("tools", {})))}
    for category in ("inputs", "scripts"):
        identity[category] = {}
        for name, raw_path in sorted(spec.get(category, {}).items()):
            path = Path(raw_path).resolve(strict=True)
            before = path.stat()
            digest = sha256_file(path)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                    after.st_size, after.st_mtime_ns, after.st_ctime_ns):
                raise ValueError("input changed while hashing: " + str(path))
            identity[category][name] = {"sha256": digest, "bytes": after.st_size}
    return identity


def key(identity):
    if identity.get("schema") != SCHEMA:
        raise ValueError("unsupported identity schema")
    namespace = identity.get("namespace", "")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", namespace) or namespace in (".", ".."):
        raise ValueError("invalid identity namespace")
    if "serialized_identity" in identity:
        # The producer supplies the exact identity bytes so its language's
        # Unicode/number rendering cannot differ from Python canonical JSON.
        serialized = identity["serialized_identity"]
        payload = json.loads(serialized)
        if payload.get("schema") != SCHEMA or payload.get("namespace") != namespace:
            raise ValueError("serialized identity namespace/schema mismatch")
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    return hashlib.sha256(canonical(identity)).hexdigest()


def location(root, identity):
    return Path(root) / "v1" / identity["namespace"] / key(identity)


def inventory(directory):
    products = {}
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError("cache products cannot contain symlinks: " + str(path))
        if path.is_dir():
            continue
        if not path.is_file():
            raise ValueError("cache products must be regular files: " + str(path))
        name = path.relative_to(directory).as_posix()
        if name == MANIFEST:
            continue
        products[name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    return products


def verify(root, identity):
    destination = location(root, identity)
    with open(destination / MANIFEST, encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("identity") != identity or manifest.get("key") != key(identity):
        raise ValueError("cache manifest identity mismatch: " + str(destination))
    if manifest.get("schema") != SCHEMA or manifest.get("state") != "complete":
        raise ValueError("cache manifest is not complete: " + str(destination))
    if not manifest.get("products") or inventory(destination) != manifest["products"]:
        raise ValueError("cache product integrity mismatch: " + str(destination))
    return destination


def fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def link_or_copy(source, target):
    try:
        os.link(source, target)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        shutil.copy2(source, target)
    return str(target)


def publish(root, identity, source):
    """Publish a complete directory by rename, serialized for the identity.

    Source products must already be closed. A concurrent successful publisher
    wins; an existing corrupt cache entry raises rather than being overwritten.
    """
    source = Path(source).resolve(strict=True)
    destination = location(root, identity)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination.parent / (key(identity) + ".lock"), "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.exists():
            return verify(root, identity)
        if (source / MANIFEST).exists():
            raise ValueError("source contains reserved manifest filename")
        source_products = inventory(source)
        if not source_products:
            raise ValueError("refusing to publish an empty cache entry")
        stage = Path(tempfile.mkdtemp(prefix="." + key(identity) + ".", dir=destination.parent))
        try:
            for name in source_products:
                target = stage / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / name, target)
                with open(target, "rb") as handle:
                    os.fsync(handle.fileno())
            products = inventory(stage)
            if products != source_products:
                raise ValueError("source products changed during publication")
            manifest = {"schema": SCHEMA, "state": "complete", "key": key(identity),
                        "identity": identity, "products": products}
            with open(stage / MANIFEST, "wb") as handle:
                handle.write(canonical(manifest) + b"\n")
                handle.flush()
                os.fsync(handle.fileno())
            for directory in sorted((p for p in stage.rglob("*") if p.is_dir()),
                                    key=lambda p: len(p.parts), reverse=True):
                fsync_directory(directory)
            fsync_directory(stage)
            os.rename(stage, destination)
            fsync_directory(destination.parent)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return destination


def run(root, identity, products, command):
    """Restore or build declared outputs; publish only after command success."""
    for name in products:
        if Path(name).is_absolute() or ".." in Path(name).parts or name == MANIFEST:
            raise ValueError("product must be a relative path inside the task directory")
    destination = location(root, identity)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Serialize misses as well as publication: two simultaneous workflows must
    # not both perform an expensive reference/coverage preparation.
    with open(destination.parent / (key(identity) + ".build.lock"), "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return run_locked(root, identity, products, command)


def run_locked(root, identity, products, command):
    destination = location(root, identity)
    if destination.exists():
        destination = verify(root, identity)
        for name in products:
            source, target = destination / name, Path(name)
            if target.exists():
                raise ValueError("cache restore target already exists: " + name)
            if source.is_dir():
                shutil.copytree(source, target, copy_function=link_or_copy)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                link_or_copy(source, target)
        return destination
    subprocess.run(command, check=True)
    # Include only explicitly declared products, never task logs or intermediates.
    with tempfile.TemporaryDirectory(prefix=".cache-products-", dir=".") as temporary:
        bundle = Path(temporary)
        for name in products:
            source, target = Path(name), bundle / name
            if source.is_dir():
                shutil.copytree(source, target, copy_function=os.link)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                os.link(source, target)
        return publish(root, identity, bundle)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    identify_parser = commands.add_parser("identify")
    identify_parser.add_argument("--spec", required=True)
    identify_parser.add_argument("--output", required=True)
    for command in ("verify", "publish", "run"):
        subparser = commands.add_parser(command)
        subparser.add_argument("--identity", required=True)
        subparser.add_argument("--root", required=True)
        if command == "publish":
            subparser.add_argument("--source", required=True)
        if command == "run":
            subparser.add_argument("--product", action="append", required=True)
            subparser.add_argument("build_command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command == "identify":
        with open(args.spec, encoding="utf-8") as handle:
            identity = identify(json.load(handle))
        Path(args.output).write_bytes(canonical(identity) + b"\n")
        print(key(identity))
    else:
        with open(args.identity, encoding="utf-8") as handle:
            identity = json.load(handle)
        if args.command == "run":
            command = args.build_command
            if command and command[0] == "--":
                command = command[1:]
            if not command:
                parser.error("run requires a build command after --")
            print(run(args.root, identity, args.product, command))
        elif args.command == "verify":
            print(verify(args.root, identity))
        else:
            print(publish(args.root, identity, args.source))


if __name__ == "__main__":
    main()
