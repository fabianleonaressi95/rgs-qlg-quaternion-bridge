import os
import yaml
import subprocess

TOTAL_NODES = 64

def generate_compose_file():
    services = {}
    for i in range(TOTAL_NODES):
        service_name = f"node_{i:02d}"
base_delay = 20 + (i % 7) * 15  # Latenza tra 20ms e 110ms
        services[service_name] = {
            "build": ".",
            "container_name": service_name,
"cap_add": ["NET_ADMIN"],
            "environment": {
                "NODE_ID": i,
                "TOTAL_NODES": TOTAL_NODES
            },
            "networks": ["aperiodic_mesh"]
        }

    compose_config = {
       
        "services": services,
        "networks": {
            "aperiodic_mesh": {
                "driver": "bridge"
            }
        }
    }

    with open("docker-compose.yml", "w") as f:
        yaml.dump(compose_config, f, sort_keys=False)
    
    print(f"[Orchestrator] Generato file docker-compose.yml con {TOTAL_NODES} nodi aperiodici.")

def run_benchmark():
    print("[Orchestrator] Avvio container e raccolta metriche di benchmark...")
    result = subprocess.run(["docker", "compose", "up", "--build", "--abort-on-container-exit"], capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("Log di sistema/errori:\n", result.stderr)

    print("[Orchestrator] Pulizia risorse Docker...")
    subprocess.run(["docker", "compose", "down", "--volumes"], stdout=subprocess.DEVNULL)
    print("[Orchestrator] Test completato.")

if __name__ == "__main__":
    generate_compose_file()
    run_benchmark()
