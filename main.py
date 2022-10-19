import subprocess
import os


class Logger(object):
    def __init__(self):
        self._log = ''

    def write(self, s, end='\n'):
        self._log += s + end

    def execute(self, command):
        run_result = subprocess.run(command, shell=True, check=False,
                                    capture_output=True)
        self._log += 'running command: "' + command + '"\n'
        self._log += 'command output: \n'
        self._log += 'stdout: \n' + ('None' if not run_result.stdout else str(run_result.stdout))
        self._log += 'stderr: \n' + ('None' if not run_result.stderr else str(run_result.stderr))

    def log(self):
        return self._log
        
        
        
        
    

raw_data_dir = './raw_data'

def get_raw_data(output_dir, logger):
    logger.write('Getting raw data...')
    if not os.path.exists(output_dir):
        logger.execute('mkdir ' + output_dir)

    with open('files_to_download.csv') as inp:
        _ = next(inp)
        for line in inp:
            download_link, filename, description = line.strip().split(',')
            if not os.path.exists(os.path.join(output_dir, filename)):
                logger.write('Getting ' + filename)
                logger.execute('wget -O {} {}'.format(os.path.join(output_dir, filename),
                                                      download_link))

    logger.write('Done.')


def main():
    logger = Logger()
    get_raw_data(raw_data_dir, logger)
    with open('log.txt') as out:
        print(logger.log, file=out)

if __name__ == "__main__":
    main()
