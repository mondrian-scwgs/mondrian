import sys
import argparse



def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id',
                        help='sample id.',
                        required=True
                        )
    parser.add_argument('--sample-ids',
                        help='List of sample ids.',
                        nargs='*',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

def main():
    args = get_args()
    for i, sample_id in enumerate(args['sample_ids']):
        if sample_id == args['sample_id']:
            sys.stdout.write(str(i))

if __name__ == "__main__":
    main() 
