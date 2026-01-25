"""
High-level, testable pipelines used by notebooks and user-facing scripts.

These functions are intentionally thin wrappers around the core `windtunnel.*`
classes and functions. They keep IO/plotting optional so future GUI work can
reuse the same compute layer.
"""


