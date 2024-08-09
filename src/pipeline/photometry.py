import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
from vso import phot
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--object', type=str, required=True, help='Object to be worked')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-c', '--common-dir', type=str, required=True, help='File tree root')
    parser.add_argument('-i', '--image-dir', type=str, default=None, help='Image directory')

    return parser.parse_args()

def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object

    session_dir = args.common_dir /  Path('session') / object_dir
    calibr_dir =  args.common_dir / 'calibr'

    source_root = args.image_dir  / object_dir / 'Light'

    p = phot.BulkPhotometry(session_dir, calibr_dir)
    p.process(source_root)

    return 0

# Example: python3 photometry.py -O RR_Lyr -t 20230704 -c /srv/public

if __name__ == '__main__':
    sys.exit(main())