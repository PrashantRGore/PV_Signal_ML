from huggingface_hub import HfApi, login
import os

print("=== Hugging Face Dataset Upload ===\n")

# Login to Hugging Face
print("Step 1: Login to Hugging Face")
print("You will be prompted for your Hugging Face token.")
print("Get your token from: https://huggingface.co/settings/tokens\n")

try:
    login()
    print("✓ Login successful!\n")
except Exception as e:
    print(f"✗ Login failed: {e}")
    print("\nCreate a token at: https://huggingface.co/settings/tokens")
    print("Then run: huggingface-cli login")
    exit(1)

# Initialize API
api = HfApi()

# Check if repository exists, if not create it
repo_id = os.getenv("HF_REPO_ID", "PrashantRGore/synthetic-faers-1m-v3")
print(f"Step 2: Creating/checking repository: {repo_id}")

try:
    api.create_repo(
        repo_id=repo_id,
        repo_type="dataset",
        exist_ok=True,
        private=False
    )
    print(f"✓ Repository ready: https://huggingface.co/datasets/{repo_id}\n")
except Exception as e:
    print(f"⚠ Note: {e}\n")

# Upload the main CSV file (PV Signal ML compatible version)
print("Step 3: Uploading synthetic_faers_1m_pv_signal_ml.csv (~384 MB)")
print("This may take 2-5 minutes depending on your internet speed...")

try:
    api.upload_file(
        path_or_fileobj="synthetic_faers_1m_pv_signal_ml.csv",
        path_in_repo="synthetic_faers_1m_pv_signal_ml.csv",
        repo_id=repo_id,
        repo_type="dataset",
    )
    print("✓ Main dataset uploaded!\n")
except Exception as e:
    print(f"✗ Upload failed: {e}\n")
    exit(1)

# Upload the original version too (optional)
print("Step 4: Uploading synthetic_faers_1m_v3.csv (original column names)")

try:
    api.upload_file(
        path_or_fileobj="synthetic_faers_1m_v3.csv",
        path_in_repo="synthetic_faers_1m_v3.csv",
        repo_id=repo_id,
        repo_type="dataset",
    )
    print("✓ Original dataset uploaded!\n")
except Exception as e:
    print(f"⚠ Original version upload failed: {e}\n")

# Upload README if it exists
print("Step 5: Uploading README.md")

if os.path.exists("README.md"):
    try:
        api.upload_file(
            path_or_fileobj="README.md",
            path_in_repo="README.md",
            repo_id=repo_id,
            repo_type="dataset",
        )
        print("✓ README uploaded!\n")
    except Exception as e:
        print(f"⚠ README upload skipped: {e}\n")
else:
    print("⚠ README.md not found - skipped\n")

print("="*60)
print("✓✓✓ UPLOAD COMPLETE ✓✓✓")
print("="*60)
print(f"\nYour dataset is now available at:")
print(f"https://huggingface.co/datasets/{repo_id}")
print(f"\nYou can load it with:")
print(f"from datasets import load_dataset")
print(f"dataset = load_dataset('{repo_id}', split='train')")
print("="*60)
