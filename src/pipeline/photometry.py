import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
from vso import phot
from vso import util
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

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    img_layout = util.ImageLayout(args.image_dir)

    p = phot.BulkPhotometry(work_layout.get_session(session).root_dir,
                            work_layout.calibr_dir)
    p.process(img_layout.get_images(session).lights_dir)

    return 0

# Example: python3 photometry.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())