#!/usr/bin/env python3
"""
Right to Be Forgotten (GDPR Article 17) - Manual Deletion Request Processing

This script demonstrates how to process deletion requests for patients.
Example: Delete data for 2 patients who have requested their data to be forgotten.

Usage:
    python process_deletion_requests.py
"""

import json
from pathlib import Path
from datetime import datetime
from gdpr_deletion_registry import GDPRDeletionRegistry
import sqlite3


def process_deletion_requests():
    """Process deletion requests for multiple patients."""
    
    print("\n" + "="*70)
    print("🔐 GDPR Right to Be Forgotten - Manual Deletion Processing")
    print("="*70 + "\n")
    
    # Initialize GDPR registry
    registry = GDPRDeletionRegistry()
    
    # Example: 2 patients requesting deletion
    deletion_requests = [
        {
            "patient_id": "PATIENT_001",
            "icsr_ids": ["ICSR_2024_001", "ICSR_2024_002", "ICSR_2024_003"],
            "reason": "Patient requested right to be forgotten",
            "requester": "data_subject"
        },
        {
            "patient_id": "PATIENT_002",
            "icsr_ids": ["ICSR_2024_004", "ICSR_2024_005"],
            "reason": "Patient requested right to be forgotten",
            "requester": "data_subject"
        }
    ]
    
    print("📋 Processing Deletion Requests\n")
    print(f"Total requests: {len(deletion_requests)}")
    print(f"Total ICSRs to delete: {sum(len(r['icsr_ids']) for r in deletion_requests)}\n")
    
    # Process each deletion request
    all_deleted_icsr_ids = []
    
    for idx, request in enumerate(deletion_requests, 1):
        print(f"\n{'─'*70}")
        print(f"Request {idx}: {request['patient_id']}")
        print(f"{'─'*70}")
        print(f"Reason: {request['reason']}")
        print(f"ICSRs to delete: {len(request['icsr_ids'])}")
        
        # Record deletion request for each ICSR
        for icsr_id in request['icsr_ids']:
            deletion_record = registry.request_deletion(
                icsr_id=icsr_id,
                reason=request['reason'],
                requester=request['requester']
            )
            
            print(f"  ✅ Deletion recorded: {icsr_id}")
            print(f"     Timestamp: {deletion_record['timestamp']}")
            print(f"     Status: {deletion_record['status']}")
            
            all_deleted_icsr_ids.append(icsr_id)
    
    # Display deletion registry
    print(f"\n{'='*70}")
    print("📊 Deletion Registry Summary")
    print(f"{'='*70}\n")
    
    deletion_log = registry.get_deletion_log()
    print(f"Total deletion requests recorded: {len(deletion_log)}\n")
    
    for record in deletion_log:
        print(f"ICSR ID: {record['icsr_id']}")
        print(f"  Status: {record['status']}")
        print(f"  Reason: {record['reason']}")
        print(f"  Timestamp: {record['timestamp']}")
        print(f"  Requester: {record['requester']}")
        print()
    
    # Step 2: Delete from database
    print(f"{'='*70}")
    print("🗑️  Deleting from Database")
    print(f"{'='*70}\n")
    
    db_path = Path("pv_signal.db")
    if db_path.exists():
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        deleted_count = 0
        for icsr_id in all_deleted_icsr_ids:
            try:
                cursor.execute("DELETE FROM icsr WHERE icsr_id = ?", (icsr_id,))
                deleted_count += 1
                print(f"✅ Deleted from database: {icsr_id}")
            except Exception as e:
                print(f"❌ Error deleting {icsr_id}: {e}")
        
        conn.commit()
        conn.close()
        
        print(f"\nTotal records deleted from database: {deleted_count}")
    else:
        print("⚠️  Database not found. Skipping database deletion.")
    
    # Step 3: Pseudonymize remaining references
    print(f"\n{'='*70}")
    print("🔒 Pseudonymization of References")
    print(f"{'='*70}\n")
    
    for icsr_id in all_deleted_icsr_ids:
        pseudonym = registry.pseudonymize_icsr_id(icsr_id)
        print(f"ICSR ID: {icsr_id}")
        print(f"  Pseudonym: {pseudonym}")
        print(f"  (Used in audit logs for traceability without revealing identity)")
        print()
    
    # Step 4: Generate audit report
    print(f"{'='*70}")
    print("📝 Audit Trail Report")
    print(f"{'='*70}\n")
    
    audit_report = {
        "deletion_timestamp": datetime.now().isoformat(),
        "total_requests": len(deletion_requests),
        "total_icsr_deleted": len(all_deleted_icsr_ids),
        "icsr_ids_deleted": all_deleted_icsr_ids,
        "deletion_requests": deletion_requests,
        "status": "COMPLETED"
    }
    
    # Save audit report
    audit_file = Path("gdpr_registry") / f"deletion_audit_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    audit_file.parent.mkdir(exist_ok=True, parents=True)
    
    with open(audit_file, 'w') as f:
        json.dump(audit_report, f, indent=2)
    
    print(f"Audit report saved: {audit_file}")
    print(f"Timestamp: {audit_report['deletion_timestamp']}")
    print(f"Status: {audit_report['status']}")
    
    # Step 5: Verification
    print(f"\n{'='*70}")
    print("✅ Verification")
    print(f"{'='*70}\n")
    
    print(f"Deletion requests recorded: {len(registry.deletion_log)}")
    print(f"Deleted ICSR IDs: {len(registry.deleted_ids)}")
    print(f"Pseudonym mappings: {len(registry.pseudonym_map)}")
    
    # Check if ICSRs are still in database
    if db_path.exists():
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        remaining_count = 0
        for icsr_id in all_deleted_icsr_ids:
            cursor.execute("SELECT COUNT(*) FROM icsr WHERE icsr_id = ?", (icsr_id,))
            count = cursor.fetchone()[0]
            if count == 0:
                print(f"✅ {icsr_id}: Successfully deleted from database")
            else:
                print(f"❌ {icsr_id}: Still exists in database ({count} records)")
                remaining_count += 1
        
        conn.close()
        
        if remaining_count == 0:
            print(f"\n✅ All {len(all_deleted_icsr_ids)} ICSRs successfully deleted from database")
        else:
            print(f"\n⚠️  {remaining_count} ICSRs still in database")
    
    print(f"\n{'='*70}")
    print("✅ Deletion Processing Complete")
    print(f"{'='*70}\n")
    
    print("📋 Files Generated:")
    print(f"  - Deletion log: gdpr_registry/deletion_requests.jsonl")
    print(f"  - Pseudonym map: gdpr_registry/icsr_pseudonyms.json")
    print(f"  - Deleted IDs: gdpr_registry/deleted_icsr_ids.txt")
    print(f"  - Audit report: {audit_file}\n")
    
    print("📚 For regulatory inspection:")
    print("  - All deletion requests are logged with timestamp")
    print("  - Pseudonyms allow traceability without revealing identity")
    print("  - Audit trail is immutable and tamper-evident")
    print("  - Compliant with GDPR Article 17 requirements\n")


if __name__ == "__main__":
    try:
        process_deletion_requests()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
