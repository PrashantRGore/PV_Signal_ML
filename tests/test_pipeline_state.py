import pytest
import pandas as pd
from src.utils.pipeline_state import PipelineState


def test_pipeline_state_current_data_truthiness():
    """Verify that current_data property does not raise DataFrame truth value error."""
    ps = PipelineState()
    assert ps.current_data is None
    
    df_raw = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
    ps.set_raw_data(df_raw, source="local")
    
    # Should not raise ValueError: The truth value of a DataFrame is ambiguous
    cur = ps.current_data
    assert cur is not None
    assert len(cur) == 2
    
    df_eng = pd.DataFrame({f"col_{i}": [i, i+1] for i in range(10)})
    ps.set_engineered_data(df_eng)
    
    cur_eng = ps.current_data
    assert cur_eng is not None
    assert len(cur_eng.columns) == 10

