import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
from vso.util import WorkLayout, FrameType
from vso.reduce import MasterBuilder


def parse_args():
    parser = argparse.ArgumentParser(
        description='Traverse image directory and create master frames'
        'from all discovered calibration images')
    parser.add_argument('-w', '--work-dir', type=str,
                        required=True, help='Work directory')
    parser.add_argument('-i', '--image-dir', type=str,
                        default=None, help='Image directory')
    parser.add_argument('-T', '--type', type=FrameType,
                        required=True, help='Calibration frame type ', choices=list(FrameType))
    parser.add_argument('-m', '--memory-limit', type=int,
                        default=None, help='memory limit in MB')
    parser.add_argument('--no-cleanup', action='store_true',
                        default=False, help='Do not remove temporary files')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()


def main():
    args = parse_args()
    layout = WorkLayout(args.work_dir)
    builder = MasterBuilder(args.type,
                            args.image_dir,
                            layout.calibr_dir,
                            layout.tmp_dir,
                            overwrite=args.overwrite,
                            delete_tmp=not args.no_cleanup)
    builder.process()


# Example: python3 create_master.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img/20240101/Calibr

if __name__ == '__main__':
    sys.exit(main())
