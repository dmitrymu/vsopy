import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
from vso import phot
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from images for the session'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-i', '--image-dir', type=str, default=None, help='Image directory')

    return parser.parse_args()

def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object

    session_dir = args.work_dir /  Path('session') / object_dir
    calibr_dir =  args.work_dir / 'calibr'

    source_root = args.image_dir  / object_dir / 'Light'

    p = phot.BulkPhotometry(session_dir, calibr_dir)
    p.process(source_root)

    return 0

# Example: python3 photometry.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())