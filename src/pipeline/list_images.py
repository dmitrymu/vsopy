import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Traverse session image directory abd create image list'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-i', '--image-dir', type=str, required=True, default=None, help='Image directory')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    img_layout = util.ImageLayout(args.image_dir)

    images = util.session_image_list(img_layout.get_images(session))
    images.meta.update({'object': args.object})

    images.write(work_layout.get_session(session).images_file_path,
                format='ascii.ecsv', overwrite=args.overwrite)

    return 0

# Example: python3 list_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())