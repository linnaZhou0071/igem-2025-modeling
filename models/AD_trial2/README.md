# AD_trial2

2D occurrence with EPS/biofilm capture. Source files are retained from the original 2025 iGEM GitLab project.

- Active configuration: [config/settings.xml](config/settings.xml)
- Initial conditions: [config/cells.csv](config/cells.csv)
- Behavior grammar: [config/cell_rules.csv](config/cell_rules.csv)
- Custom dynamics: [custom_modules/custom.cpp](custom_modules/custom.cpp)

From the repository root:

```bash
python3 scripts/run_model.py AD_trial2 --threads 2 --max-time 1
```

Remove `--max-time 1` to use the full configured horizon. Read the [reproduction guide](../../docs/REPRODUCIBILITY.md) for configuration selection, historical path corrections and validation limits.
