import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Make batches from image list'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    # parser.add_argument('-i', '--image-dir', type=str, required=True, default=None, help='Image directory')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    session_layout = util.WorkLayout(args.work_dir).get_session(session)

    batches, batch_images = util.batch_session_images(session_layout.images_file_path)
    batches.meta.update({'object': args.object})

    batches.write(session_layout.batches_file_path,
                  format='ascii.ecsv', overwrite=args.overwrite)

    batch_images.write(session_layout.batch_images_file_path,
                       format='ascii.ecsv', overwrite=args.overwrite)

    return 0

# Example: python3 batch_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())