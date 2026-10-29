import pandas as pd
import numpy as np
import random
from datetime import datetime

np.random.seed(42)
random.seed(42)

# === STEP 1: Define 10 signal drug profiles ===
SIGNAL_DRUGS = {
    'ANTICOAG-XR': {
        'total_cases': 40_000,
        'signal_events': {
            'Gastrointestinal haemorrhage': 0.09,
            'Cerebral haemorrhage': 0.03,
            'Haematuria': 0.04
        },
        'common_aes': ['Nausea', 'Headache', 'Dizziness', 'Fatigue', 'Abdominal pain'],
        'common_ae_rate': 0.50,
        'background_rate': 0.34
    },
    
    'Hepatosan': {
        'total_cases': 35_000,
        'signal_events': {
            'Hepatic failure': 0.08,
            'Hepatotoxicity': 0.06,
            'Jaundice': 0.05
        },
        'common_aes': ['Nausea', 'Abdominal pain', 'Fatigue', 'Vomiting', 'Anorexia'],
        'common_ae_rate': 0.52,
        'background_rate': 0.29
    },
    
    'NEPHRO-GUARD': {
        'total_cases': 38_000,
        'signal_events': {
            'Acute kidney injury': 0.10,
            'Renal failure': 0.07,
            'Haematuria': 0.05
        },
        'common_aes': ['Nausea', 'Headache', 'Oedema peripheral', 'Fatigue', 'Dizziness'],
        'common_ae_rate': 0.48,
        'background_rate': 0.30
    },
    
    'CARDIO-STABIL': {
        'total_cases': 42_000,
        'signal_events': {
            'Cardiac failure': 0.08,
            'Arrhythmia': 0.06,
            'Myocardial infarction': 0.04
        },
        'common_aes': ['Chest pain', 'Dyspnoea', 'Fatigue', 'Dizziness', 'Palpitations'],
        'common_ae_rate': 0.55,
        'background_rate': 0.27
    },
    
    'NEURO-CALM': {
        'total_cases': 33_000,
        'signal_events': {
            'Seizure': 0.09,
            'Tremor': 0.07,
            'Confusion': 0.05
        },
        'common_aes': ['Headache', 'Dizziness', 'Drowsiness', 'Nausea', 'Insomnia'],
        'common_ae_rate': 0.50,
        'background_rate': 0.29
    },
    
    'ONCO-SHIELD': {
        'total_cases': 30_000,
        'signal_events': {
            'Bone marrow failure': 0.10,
            'Febrile neutropenia': 0.08,
            'Thrombocytopenia': 0.06
        },
        'common_aes': ['Nausea', 'Vomiting', 'Fatigue', 'Alopecia', 'Anaemia'],
        'common_ae_rate': 0.60,
        'background_rate': 0.16
    },
    
    'IMMUNO-BOOST': {
        'total_cases': 28_000,
        'signal_events': {
            'Anaphylactic reaction': 0.07,
            'Hypersensitivity': 0.06,
            'Angioedema': 0.05
        },
        'common_aes': ['Pyrexia', 'Rash', 'Pruritus', 'Headache', 'Nausea'],
        'common_ae_rate': 0.53,
        'background_rate': 0.29
    },
    
    'PSYCH-EASE': {
        'total_cases': 32_000,
        'signal_events': {
            'Suicidal ideation': 0.08,
            'Aggression': 0.05,
            'Hallucination': 0.04
        },
        'common_aes': ['Insomnia', 'Anxiety', 'Headache', 'Nausea', 'Dizziness'],
        'common_ae_rate': 0.51,
        'background_rate': 0.32
    },
    
    'DIABETES-CTRL': {
        'total_cases': 45_000,
        'signal_events': {
            'Hypoglycaemia': 0.12,
            'Diabetic ketoacidosis': 0.05,
            'Hyperglycaemia': 0.04
        },
        'common_aes': ['Nausea', 'Headache', 'Dizziness', 'Abdominal pain', 'Diarrhoea'],
        'common_ae_rate': 0.54,
        'background_rate': 0.25
    },
    
    'RESPIRO-CLEAR': {
        'total_cases': 36_000,
        'signal_events': {
            'Pneumonia': 0.09,
            'Bronchospasm': 0.06,
            'Respiratory failure': 0.04
        },
        'common_aes': ['Cough', 'Dyspnoea', 'Headache', 'Nausea', 'Throat irritation'],
        'common_ae_rate': 0.52,
        'background_rate': 0.29
    }
}

# Background event pool (creates diversity for robust PRR)
BACKGROUND_EVENTS = [
    'Pyrexia', 'Rash', 'Pruritus', 'Asthenia', 'Malaise', 'Chest pain', 'Back pain',
    'Dyspnoea', 'Cough', 'Oedema peripheral', 'Hypertension', 'Hypotension', 
    'Tachycardia', 'Bradycardia', 'Constipation', 'Diarrhoea', 'Flatulence',
    'Dyspepsia', 'Gastritis', 'Urinary tract infection', 'Dysuria', 'Arthralgia',
    'Myalgia', 'Pain in extremity', 'Muscle spasms', 'Depression', 'Anxiety',
    'Insomnia', 'Somnolence', 'Vertigo', 'Paraesthesia', 'Hypoaesthesia',
    'Vision blurred', 'Tinnitus', 'Dry mouth', 'Increased appetite', 'Decreased appetite',
    'Weight increased', 'Weight decreased', 'Hot flush', 'Chills', 'Night sweats',
    'Peripheral coldness', 'Feeling abnormal', 'Feeling jittery', 'Hunger', 'Thirst',
    'Pollakiuria', 'Polyuria', 'Nocturia', 'Micturition urgency', 'Haematochezia',
    'Melaena', 'Haematemesis', 'Rectal haemorrhage', 'Epistaxis', 'Gingival bleeding',
    'Petechiae', 'Purpura', 'Ecchymosis', 'Contusion', 'Swelling face',
    'Swelling', 'Joint swelling', 'Oedema', 'Fluid retention', 'Ascites',
    'Pleural effusion', 'Photosensitivity reaction', 'Skin discolouration', 'Alopecia',
    'Hirsutism', 'Nail disorder', 'Skin disorder', 'Urticaria', 'Dermatitis',
    'Eczema', 'Erythema', 'Blister', 'Acne', 'Seborrhoea', 'Hyperhidrosis'
]

# 15 filler drugs (no known signals)
FILLER_DRUGS = [f"DRUG-{i:02d}" for i in range(1, 16)]

# === STEP 2: Generate cases ===
def generate_v4_quick(target_size=1_000_000):
    print(f"🚀 Generating {target_size:,} cases for PV_Signal_ML v4...")
    print(f"⏰ Started at: {datetime.now().strftime('%H:%M:%S')}")
    
    records = []
    case_id = 1
    
    # Allocate 70% to signal drugs, 30% to filler drugs
    signal_allocation = int(target_size * 0.70)
    filler_allocation = int(target_size * 0.30)
    
    # === Generate signal drug cases ===
    total_signal_cases = sum(d['total_cases'] for d in SIGNAL_DRUGS.values())
    
    for drug_name, profile in SIGNAL_DRUGS.items():
        n_cases = int(profile['total_cases'] * (signal_allocation / total_signal_cases))
        print(f"  📊 {drug_name}: {n_cases:,} cases")
        
        for _ in range(n_cases):
            # Choose event type
            rand = random.random()
            
            if rand < sum(profile['signal_events'].values()):
                # Signal event
                event = random.choices(
                    list(profile['signal_events'].keys()),
                    weights=list(profile['signal_events'].values())
                )[0]
                outcome = random.choice(['Hospitalization', 'Life-threatening', 'Death'])
                
            elif rand < sum(profile['signal_events'].values()) + profile['common_ae_rate']:
                # Common AE
                event = random.choice(profile['common_aes'])
                outcome = random.choice(['Other', 'Other', 'Hospitalization'])
                
            else:
                # Background event (creates diversity for robust PRR)
                event = random.choice(BACKGROUND_EVENTS)
                outcome = 'Other'
            
            # Generate minimal ICSR
            record = {
                'case_id': f"US-{case_id:09d}",
                'suspect_drug': drug_name,
                'event_pt': event,
                'age': random.randint(18, 85),
                'gender': random.choice(['M', 'F']),
                'reporter': random.choice(['Physician', 'Pharmacist', 'Consumer']),
                'outcome': outcome,
                'report_year': random.randint(2019, 2024),
                'country': random.choice(['US', 'GB', 'DE', 'FR', 'CA'])
            }
            records.append(record)
            case_id += 1
    
    # === Generate filler drug cases (pure noise) ===
    filler_per_drug = filler_allocation // len(FILLER_DRUGS)
    
    print(f"\n  📊 Filler drugs ({len(FILLER_DRUGS)} drugs): {filler_per_drug:,} cases each")
    
    for drug_name in FILLER_DRUGS:
        for _ in range(filler_per_drug):
            event = random.choice(BACKGROUND_EVENTS)
            
            record = {
                'case_id': f"US-{case_id:09d}",
                'suspect_drug': drug_name,
                'event_pt': event,
                'age': random.randint(18, 85),
                'gender': random.choice(['M', 'F']),
                'reporter': random.choice(['Physician', 'Pharmacist', 'Consumer']),
                'outcome': random.choice(['Other', 'Hospitalization']),
                'report_year': random.randint(2019, 2024),
                'country': random.choice(['US', 'GB', 'DE', 'FR', 'CA'])
            }
            records.append(record)
            case_id += 1
    
    df = pd.DataFrame(records)
    print(f"\n✅ Generated {len(df):,} cases")
    print(f"⏰ Completed at: {datetime.now().strftime('%H:%M:%S')}")
    
    # Quick validation
    print("\n🔍 Validation Check:")
    for drug in list(SIGNAL_DRUGS.keys())[:3]:
        drug_df = df[df['suspect_drug'] == drug]
        unique_events = drug_df['event_pt'].nunique()
        status = "✅" if unique_events >= 50 else "⚠️"
        print(f"  {status} {drug}: {len(drug_df):,} cases, {unique_events} unique events")
    
    return df

# === RUN ===
if __name__ == "__main__":
    df_v4 = generate_v4_quick(target_size=1_000_000)
    
    output_path = r"C:\Users\koreo\PV_Signal_ML\data\synthetic_faers_1m_v4.csv"
    print(f"\n💾 Saving to: {output_path}")
    df_v4.to_csv(output_path, index=False)
    
    print(f"\n✅ DATASET GENERATION COMPLETE!")
    print(f"📊 Columns: {list(df_v4.columns)}")
    print(f"📈 Shape: {df_v4.shape}")
    print(f"💽 File size: ~{df_v4.memory_usage(deep=True).sum() / 1024**2:.1f} MB in memory")
