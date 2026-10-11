#!/usr/bin/env python3
"""Build an isolated PhysiCell project and run a portable copy of its settings."""
import argparse
import shutil
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
MODELS = ('AD_trial1', 'AD_trial2', 'Recurrence', 'Treatment', 'AD_3D', 'Recurrence_3D')
ACTIVE_SETTINGS = {
    'AD_trial1': 'settings', 'AD_trial2': 'settings.xml',
    'Recurrence': 'settings.xml', 'Treatment': 'settings.xml',
    'AD_3D': 'PhysiCell_settings.xml', 'Recurrence_3D': 'PhysiCell_settings.xml',
}


def prepare(name, threads, max_time=None):
    source = ROOT / 'models' / name
    build = ROOT / 'build' / name
    build.mkdir(parents=True, exist_ok=True)
    for directory in ('BioFVM', 'core', 'modules', 'licenses'):
        shutil.copytree(ROOT / 'vendor/PhysiCell' / directory, build / directory, dirs_exist_ok=True)
    for filename in ('VERSION.txt', 'CITATION.txt'):
        shutil.copy2(ROOT / 'vendor/PhysiCell' / filename, build / filename)
    for filename in ('main.cpp', 'Makefile'):
        shutil.copy2(source / filename, build / filename)
    for directory in ('custom_modules', 'config'):
        shutil.copytree(source / directory, build / directory, dirs_exist_ok=True)
    colors = source / 'custom_modules/cell_color_map.txt'
    if colors.exists():
        shutil.copy2(colors, build / colors.name)

    # Studio saved some active settings under a different filename.
    active = source / 'config' / ACTIVE_SETTINGS[name]
    tree = ET.parse(active)
    settings = tree.getroot()
    for parent in (settings.find('initial_conditions/cell_positions'),
                   *settings.findall('cell_rules/rulesets/ruleset')):
        if parent is not None:
            parent.find('folder').text = './config'
    settings.find('save/folder').text = './output'
    settings.find('parallel/omp_num_threads').text = str(threads)
    if max_time is not None:
        settings.find('overall/max_time').text = str(max_time)
    # Replays should start from the supplied CSV, without legacy machine paths.
    output = build / 'config/portable.xml'
    tree.write(output, encoding='utf-8', xml_declaration=True)
    (build / 'output').mkdir(exist_ok=True)
    return build, output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', choices=MODELS)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--jobs', type=int, default=2, help='Concurrent compiler processes')
    parser.add_argument('--max-time', type=float, help='Override simulated duration, in minutes')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--build-only', action='store_true')
    args = parser.parse_args()
    if args.threads < 1 or args.jobs < 1 or (args.max_time is not None and args.max_time < 0):
        parser.error('threads/jobs must be positive and max-time must be nonnegative')
    build, config = prepare(args.model, args.threads, args.max_time)
    print(f'Prepared {args.model}: {config}', flush=True)
    if args.prepare_only:
        return
    subprocess.run(['make', f'-j{args.jobs}'], cwd=build, check=True)
    if not args.build_only:
        subprocess.run([str(build / 'project'), str(config)], cwd=build, check=True)


if __name__ == '__main__':
    main()
