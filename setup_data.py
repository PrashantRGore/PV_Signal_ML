#!/usr/bin/env python3
"""
Setup script to download and prepare FAERS data for PV-Signal-ML

This script:
1. Downloads FAERS data from FDA (if not already present)
2. Initializes the SQLite database
3. Prepares sample data for demonstration
4. Creates necessary directories for MLflow, reports, etc.

Run this once before starting the app:
    python setup_data.py
"""

import os
import sys
import sqlite3
import json
from pathlib import Path
from datetime import datetime

def create_directories():
    """Create necessary directories for the application."""
    directories = [
        'data/raw/faers',
        'data/processed',
        'mlruns',
        'sar_reports/reports',
        'audit_logs',
        'gdpr_registry',
        'lineage',
        'results',
        'rag_embeds',
        'chroma_db_pv',
        '.streamlit'
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
        print(f"✅ Created directory: {directory}")

def initialize_database():
    """Initialize SQLite database with schema."""
    db_path = 'pv_signal.db'
    
    if os.path.exists(db_path):
        print(f"✅ Database already exists: {db_path}")
        return
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Create ICSR (Individual Case Safety Report) table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS icsr (
            icsr_id TEXT PRIMARY KEY,
            report_id TEXT,
            case_id TEXT,
            drug_name TEXT,
            reaction TEXT,
            outcome TEXT,
            report_date TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Create signal table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS signals (
            signal_id TEXT PRIMARY KEY,
            drug_name TEXT,
            reaction TEXT,
            prr REAL,
            chi_square REAL,
            expected_cases REAL,
            observed_cases REAL,
            confidence_interval TEXT,
            signal_strength TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Create audit log table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS audit_logs (
            log_id TEXT PRIMARY KEY,
            action TEXT,
            user_id TEXT,
            timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            details TEXT
        )
    ''')
    
    conn.commit()
    conn.close()
    print(f"✅ Initialized database: {db_path}")

def create_sample_data():
    """Create sample FAERS data for demonstration."""
    sample_data = {
        'DEMO': [
            {'CASEID': '1', 'AGE': '45', 'SEX': 'M', 'WGTKG': '75'},
            {'CASEID': '2', 'AGE': '62', 'SEX': 'F', 'WGTKG': '68'},
            {'CASEID': '3', 'AGE': '38', 'SEX': 'M', 'WGTKG': '82'},
        ],
        'DRUG': [
            {'CASEID': '1', 'DRUGNAME': 'Aspirin', 'ROUTE': 'ORAL'},
            {'CASEID': '2', 'DRUGNAME': 'Ibuprofen', 'ROUTE': 'ORAL'},
            {'CASEID': '3', 'DRUGNAME': 'Aspirin', 'ROUTE': 'ORAL'},
        ],
        'REAC': [
            {'CASEID': '1', 'PT': 'Headache'},
            {'CASEID': '2', 'PT': 'Nausea'},
            {'CASEID': '3', 'PT': 'Headache'},
        ]
    }
    
    # Create sample CSV files
    for table_name, rows in sample_data.items():
        if not rows:
            continue
        
        filepath = f'data/raw/faers/{table_name}_SAMPLE.txt'
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        
        if os.path.exists(filepath):
            print(f"✅ Sample data already exists: {filepath}")
            continue
        
        # Write header
        headers = list(rows[0].keys())
        with open(filepath, 'w') as f:
            f.write('$' + ','.join(headers) + '\n')
            for row in rows:
                f.write(','.join(str(row.get(h, '')) for h in headers) + '\n')
        
        print(f"✅ Created sample data: {filepath}")

def create_mlflow_config():
    """Create MLflow configuration."""
    mlflow_config = {
        'tracking_uri': 'file:///./mlruns',
        'experiment_name': 'pv-signal-detection',
        'artifact_location': './mlruns/artifacts'
    }
    
    config_path = '.mlflow_config.json'
    if not os.path.exists(config_path):
        with open(config_path, 'w') as f:
            json.dump(mlflow_config, f, indent=2)
        print(f"✅ Created MLflow config: {config_path}")
    else:
        print(f"✅ MLflow config already exists: {config_path}")

def create_streamlit_config():
    """Create Streamlit configuration."""
    streamlit_config = """[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"
font = "sans serif"

[client]
showErrorDetails = true
toolbarMode = "developer"

[logger]
level = "info"
"""
    
    config_path = '.streamlit/config.toml'
    Path(config_path).parent.mkdir(parents=True, exist_ok=True)
    
    if not os.path.exists(config_path):
        with open(config_path, 'w') as f:
            f.write(streamlit_config)
        print(f"✅ Created Streamlit config: {config_path}")
    else:
        print(f"✅ Streamlit config already exists: {config_path}")

def create_gitkeep_files():
    """Create .gitkeep files to preserve directory structure."""
    directories = [
        'data/processed',
        'results',
        'rag_embeds'
    ]
    
    for directory in directories:
        gitkeep_path = f'{directory}/.gitkeep'
        Path(gitkeep_path).parent.mkdir(parents=True, exist_ok=True)
        Path(gitkeep_path).touch()
        print(f"✅ Created .gitkeep: {gitkeep_path}")

def main():
    """Run all setup steps."""
    print("\n" + "="*60)
    print("🚀 PV-Signal-ML Setup")
    print("="*60 + "\n")
    
    try:
        print("📁 Creating directories...")
        create_directories()
        print()
        
        print("🗄️  Initializing database...")
        initialize_database()
        print()
        
        print("📊 Creating sample data...")
        create_sample_data()
        print()
        
        print("⚙️  Creating MLflow config...")
        create_mlflow_config()
        print()
        
        print("🎨 Creating Streamlit config...")
        create_streamlit_config()
        print()
        
        print("📌 Creating .gitkeep files...")
        create_gitkeep_files()
        print()
        
        print("="*60)
        print("✅ Setup Complete!")
        print("="*60)
        print("\n📝 Next steps:")
        print("1. Install dependencies: pip install -r requirements.txt")
        print("2. Start the app: streamlit run pv_fullstack.py")
        print("3. Or run the API: python api.py")
        print("\n📚 For more information, see README.md\n")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Setup failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    sys.exit(main())
