import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
import concurrent.futures as cf
from astropy.table import QTable, vstack
from vso import phot
from vso import reduce
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from session images'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-p', '--parallel', type=int, default=1, help='Number of parallel workers')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite output files')

    return parser.parse_args()

CONTEXT = None

def make_context(args):
    print('making context')
    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)
    solver = lambda path: reduce.astap_solver(path, session_layout.solved_dir)
    matcher = reduce.CalibrationMatcher(work_layout.calibr_dir)
    centroids = QTable.read(session_layout.centroid_file_path)
    settings = util.Settings(session_layout.settings_file_path)
    global CONTEXT
    CONTEXT = (matcher, solver, centroids, settings.aperture)

def measure_image(id, path):
    print(f'measure {path}')
    global CONTEXT
    matcher, solver, centroids, aperture = CONTEXT
    result = phot.process_image(path, matcher, solver, lambda _: centroids, aperture)
    result['image_id'] = id
    return result['image_id', 'auid', 'M', 'flux', 'snr', 'peak']

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)

    images = QTable.read(session_layout.images_file_path)

    blacklist = util.Blacklist(session_layout.blacklist_file_path)

    with cf.ProcessPoolExecutor(initializer=make_context,
                                    initargs=(args,),
                                    max_workers=args.parallel) as executor:
        futures = [(image['path'], executor.submit(measure_image, image['image_id'], image['path']))
                  for image in images if not blacklist.contains(image['path'])]

    def get_result(image_result):
        path = None
        try:
            path, future = image_result
            return future.result()
        except Exception as e:
            blacklist.add(path, e)
            return None

    tables = [get_result(f) for f in futures]

    result = vstack([t for t in tables if t is not None])

    result.write(session_layout.measured_file_path, format='ascii.ecsv', overwrite=args.overwrite)
    blacklist.save(session_layout.blacklist_file_path)

    return 0

# Example: python3 measure_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())