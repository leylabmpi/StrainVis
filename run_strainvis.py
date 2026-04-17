import subprocess
import argparse
import time

RESTART_DELAY = 2  # seconds
MAX_RETRIES = 3
retry_count = 0


def run(panel_args):
    while True:
        print("\n\nStarting Panel server...")

        process = subprocess.Popen(panel_args)

        process.wait()  # wait until it exits

        exit_code = process.returncode
        print(f"\nPanel exited with code {exit_code}")

        # Optional: only restart on "intentional" exit
        if exit_code == 0:
            print("\nRestarting Panel server...")
            retry_count = 0
            time.sleep(RESTART_DELAY)
        else:
            retry_count += 1
            print(f"\nPanel crash detected (exit {exit_code}), retry {retry_count}/{MAX_RETRIES}")

            if retry_count >= MAX_RETRIES:
                print("\nToo many crashes, giving up...")
                break

            time.sleep(RESTART_DELAY)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=5005)
    args = parser.parse_args()
    app = "strain_vis.py"

    run([
        "panel", "serve", app,
        "--port", str(args.port),
        "--websocket-max-message-size", "524288000"
        ])
