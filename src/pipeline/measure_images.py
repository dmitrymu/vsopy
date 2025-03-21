import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
from astropy.table import QTable, vstack
from vso import phot
from vso import reduce
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from images for the session'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)

    solver = lambda path: reduce.astap_solver(path, session_layout.solved_dir)
    matcher = reduce.CalibrationMatcher(work_layout.calibr_dir)
    centroids = QTable.read(session_layout.centroid_file_path)
    settings = util.Settings(session_layout.settings_file_path)

    images = QTable.read(session_layout.images_file_path)

    blacklist = util.Blacklist(session_layout.blacklist_file_path)

    def measure_image(id, path):
        try:
            result = phot.process_image(path, matcher, solver, centroids, settings.aperture)
            result['image_id'] = id
            return result['image_id', 'auid', 'M', 'flux', 'snr', 'peak']
        except Exception as e:
            blacklist.add(path, e)
            return None

    per_image = [measure_image(image['image_id'], image['path'])
                 for image in images if not blacklist.contains(image['path'])]

    result = vstack([image for image in per_image if image is not None])

    result.write(session_layout.measured_file_path, format='ascii.ecsv', overwrite=args.overwrite)
    blacklist.save(session_layout.blacklist_file_path)

    return 0

# Example: python3 photometry.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())