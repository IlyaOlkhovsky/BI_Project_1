import subprocess
import os
import argparse


class Logger(object):
    log_file = 'log.txt'
    
    def __init__(self):
        if os.path.exists('log.txt'):
            os.remove('log.txt')

    def write(self, s, end='\n'):
        print(s, end=end)
        with open(self.log_file, 'a') as f_out:
            f_out.write(s + end)

    def execute(self, command, can_crash=False):
        self.write('running command: "' + command + '"\n')
        
        try:
            run_result = subprocess.run(command, shell=True, check=not can_crash,
                                        text=True, capture_output=True)
            
        except subprocess.CalledProcessError as err:
            self.write('Command exited with a nonzero return code!')
            if err.stdout:
                self.write('stdout: \n' + err.stdout)
            if err.stderr:
                self.write('stderr: \n' + err.stderr)
                
            self.write('Cannot proceed. Aborting.')
            exit(0)
            
        if run_result.stdout:
            self.write('stdout: \n' + run_result.stdout)
        if run_result.stderr:
            self.write('stderr: \n' + run_result.stderr)
            
        return (run_result.stdout, run_result.stderr)


def get_raw_data(csv_file_list, output_dir, logger):
    logger.write('-----Getting raw data...\n')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    with open(csv_file_list) as inp:
        _ = next(inp)
        for line in inp:
            download_link, filename, description = line.strip().split(',')
            logger.write('Getting ' + filename + '\n')
            
            full_filename = os.path.join(output_dir, filename)
            if os.path.exists(full_filename):
                logger.write('File already exists, passing\n')
            else:
                if not os.path.exists(full_filename + '.gz'):
                    logger.execute('wget -O {} {}'.format(full_filename + '.gz',
                                                          download_link))
                logger.execute('gzip -d {}'.format(full_filename + '.gz'))

    logger.write('Done.\n')


def inspect_raw_data(data_dir, fastq_files, logger):
    logger.write('----Inspecting raw data...')
    for f in fastq_files:
        command_output, _ = logger.execute('wc -l {}'.format(os.path.join(data_dir, f)))
        reads_number = int(command_output.split()[0]) // 4
        logger.write('Number of reads in {}: {}\n'.format(f, reads_number))

    logger.write('Done.\n')

def inspect_raw_data_fastq(data_dir, fastq_files, logger):
    logger.write('----Inspecting raw data with fastq...')
    command_output, _ = logger.execute('which fastqc', can_crash=True)
    if not command_output:
        logger.write('fastqc not found. You need to install it. Cannot proceed. Aborting')
        exit(0)

    if not os.path.exists('fastqc_output'):
        os.mkdir('fastqc_output')

    for f in fastq_files:
        if os.path.exists(os.path.join('fastqc_output', os.path.splitext(f)[0] + '_fastqc.html')):
            logger.write('File {} is already processed'.format(f))
        else:
            logger.execute('fastqc -o ./fastqc_output {}'.format(os.path.join(data_dir, f)))

    logger.write('Done.\n')
        
def trim_reads(data_dir, fastq_files, logger, trimmomatic_path):
    logger.write('-----Trimming reads with Trimmomatic...')

    # checking that necessary files exist
    command_output, _ = logger.execute('which trimmomatic', can_crash=True)
    if not command_output and not trimmomatic_path:
        logger.write('trimmomatic not found. You need to install it. Cannot proceed. Aborting')
        exit(0)


    if not os.path.exists('./TruSeq3-PE.fa') and not trimmomatic_path:
        logger.write('file TruSeq3-PE.fa not found. Aborting')
        

    # configuring trimmomatic arguments
    if os.path.exists('./TruSeq3-PE.fa'):
        adaptor_path = 'TruSeq3-PE.fa'
    else:
        adaptor_path = os.path.join(os.path.split(trimmomatic_path)[0],
                                    'adaptors', 'TruSeq3-PE.fa')

    adaptor = 'ILLUMINACLIP:{}:2:30:10'.format(adaptor_path)
    

    if trimmomatic_path:
        trimmomatic_runner = 'java -jar {}'.format(trimmomatic_path)
    else:
        trimmomatic_runner = 'trimmomatic'

        
    input_files = ' '.join(os.path.join(data_dir, input_file)
                           for input_file in fastq_files)
    
    output_files = ' '.join(os.path.join('./trimmomatic_output', output_file)
                            for output_file in ['output_forward_paired.fq.gz',
                                                'output_forward_unpaired.fq.gz',
                                                'output_reverse_paired.fq.gz',
                                                'output_reverse_unpaired.fq.gz'])
        
    format_string = ' '.join(['{trimmomatic_runner} PE -phred33',
                              '{input_files}',
                              '{output_files}',
                              '{adaptor}',
                              'LEADING:{leading}',
                              'TRAILING:{trailing}',
                              'SLIDINGWINDOW:{sliding_window_1}:{sliding_window_2}',
                              'MINLEN:{minlen}'])

    logger.execute(format_string.format(trimmomatic_runner=trimmomatic_runner,
                                        input_files=input_files,
                                        output_files=output_files,
                                        adaptor=adaptor,
                                        leading=20,
                                        trailing=20,
                                        sliding_window_1=10,
                                        sliding_window_2=20,
                                        minlen=20))
                                                                                                                                                                                                    
def get_command_line_args():
    parser = argparse.ArgumentParser(description="Script for finding mutations in DNA")
    parser.add_argument('--trimmomatic', default=None,
                        help='path to the trimmomatic jar file. If not specified, "trimmomatic" command is used')
    args =  parser.parse_args()
    return args


def main():
    args = get_command_line_args()
    raw_data_dir = './raw_data'
    
    if not os.path.exists('files_to_download.csv'):
        print('Cannot find file "files_to_download.csv". Cannot Proceed. Aborting')
        
    with open('files_to_download.csv') as f:
        lines_words = iter(line.strip().split(',') for line in f)
        next(lines_words)
        fastq_files = [filename for _, filename, _ in lines_words
                       if os.path.splitext(filename)[1] == '.fastq']
        
    logger = Logger()
    get_raw_data('files_to_download.csv', raw_data_dir, logger)

    inspect_raw_data(raw_data_dir, fastq_files, logger)
    inspect_raw_data_fastq(raw_data_dir, fastq_files, logger)

    trim_reads(raw_data_dir, fastq_files, logger, args.trimmomatic)
    inspect_trimmed_data_fastq(raw_data_dir, fastq_files, logger)
    

if __name__ == "__main__":
    main()
