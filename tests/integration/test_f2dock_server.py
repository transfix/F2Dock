#!/usr/bin/env python3

import argparse
import pathlib
import shutil
import socket
import subprocess
import sys
import time
import xmlrpc.client

SUBMITTED = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Integration smoke test for F2DockServer")
    parser.add_argument("--server-exe", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--platform", type=int, default=0)
    parser.add_argument("--startup-timeout", type=float, default=15.0)
    return parser.parse_args()


def allocate_port() -> int:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port


def terminate_process(proc: subprocess.Popen[str]) -> None:
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=5)


def fail(message: str, proc: subprocess.Popen[str] | None = None) -> None:
    details = [message]
    if proc is not None and proc.poll() is not None:
        details.append(f"server exit code: {proc.returncode}")
    raise SystemExit("\n\n".join(details))


def main() -> int:
    args = parse_args()
    server_exe = pathlib.Path(args.server_exe)
    if not server_exe.exists():
        raise SystemExit(f"F2DockServer executable not found: {server_exe}")

    work_dir = pathlib.Path(args.work_dir)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    port = allocate_port()
    proc = subprocess.Popen(
        [str(server_exe), str(port), str(args.platform)],
        cwd=work_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        text=True,
    )

    old_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(5.0)
    proxy = xmlrpc.client.ServerProxy(f"http://127.0.0.1:{port}/RPC2", allow_none=True)
    deadline = time.time() + args.startup_timeout
    methods = None
    last_error = None
    try:
        while time.time() < deadline:
            if proc.poll() is not None:
                fail(f"F2DockServer exited early with code {proc.returncode}", proc)
            try:
                methods = proxy.system.listMethods()
                break
            except Exception as exc:  # pragma: no cover - exercised only during startup retries
                last_error = exc
                time.sleep(0.25)

        if methods is None:
            fail(f"F2DockServer did not start listening in time: {last_error}", proc)

        required_methods = {
            "Submit",
            "GetOutput",
            "getJobStatus",
            "system.listMethods",
            "system.methodHelp",
            "system.multicall",
        }
        missing = sorted(required_methods.difference(methods))
        if missing:
            fail(f"F2DockServer missing expected XML-RPC methods: {missing}", proc)

        platform_ini = work_dir / "platform.ini"
        if not platform_ini.exists():
            fail("F2DockServer did not create platform.ini", proc)
        platform_value = platform_ini.read_text(encoding="utf-8").strip()
        if platform_value != str(args.platform):
            fail(f"Unexpected platform.ini contents: {platform_value!r}", proc)

        status = proxy.getJobStatus("99999")
        if status != SUBMITTED:
            fail(f"Expected getJobStatus on unknown job to return {SUBMITTED}, got {status}", proc)

        method_help = proxy.system.methodHelp("Submit")
        if not isinstance(method_help, str):
            fail("system.methodHelp('Submit') did not return a string", proc)

    finally:
        socket.setdefaulttimeout(old_timeout)
        terminate_process(proc)

    return 0


if __name__ == "__main__":
    sys.exit(main())
