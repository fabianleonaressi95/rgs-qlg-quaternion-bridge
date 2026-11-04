import os

def unpack_mu_file(mu_path="mu", output_dir="."):
    if not os.path.exists(mu_path):
        print(f"Errore: Il file '{mu_path}' non è stato trovato.")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    with open(mu_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Modifica il delimitatore qui sotto in base a come hai unito i file in Kaggle
    # Esempio tipico: === FILE: nome_script.py ===
    delimiter = "=== FILE: "
    chunks = content.split(delimiter)
    
    if len(chunks) <= 1:
        print("Nessun delimitatore trovato con '=== FILE: '. Controlla l'intestazione del file mu.")
        return

    for chunk in chunks[1:]:
        lines = chunk.splitlines()
        filename = lines[0].replace(" ===", "").strip()
        file_content = "\n".join(lines[1:])
        
        out_path = os.path.join(output_dir, filename)
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        
        with open(out_path, "w", encoding="utf-8") as out_f:
            out_f.write(file_content.strip() + "\n")
        print(f"[RIPRISTINATO] -> {out_path}")

if __name__ == "__main__":
    unpack_mu_file()
