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


def get_raw_data(file_list, output_dir, logger):
    logger.write('-----Getting raw data...\n')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    for filename in file_list:
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

def inspect_data_fastq(data_dir, fastq_files, logger, log_message):
    logger.write(log_message)
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
        logger.write('\n'.join(['trimmomatic not found.',
                                'You need to install it or specify the path to a jar file.',
                                'Cannot proceed. Aborting']))
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
                            for output_file in ['output_forward_paired.fq',
                                                'output_forward_unpaired.fq',
                                                'output_reverse_paired.fq',
                                                'output_reverse_unpaired.fq'])
    if all(map(os.path.exists, output_files.split())):
        logger.write('Output files already exist. Passing.')
    else:
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
    logger.write('Done.')
                                                                                                                                                                                                    
def check_trimmed_reads_count(trimmomatic_output_dir, logger):
    logger.write('-----Checking trimmed reads count...')
    logger.execute('zcat {} | wc -l'.format(os.path.join(trimmomatic_output_dir,
                                                         'output_forward_paired.fq.gz')))

def index_reference_file(ref_genome_file, logger):
    logger.write('-----Indexing reference file...')
    command_output, _ = logger.execute('which bwa', can_crash=True)
    if not command_output:
        logger.write('bwa not found. You need to install it. Cannot proceed. Aborting')
        exit(0)
    
    logger.execute('bwa index {}'.format(ref_genome_file))
    logger.write('Done.')



def align_reads(ref_genome_file, reads_files, output_file, logger):
    logger.write('-----Aligning reads...')
    if os.path.exists(output_file):
        logger.write('File {} already exists. Passing'.format(output_file))
    else:
        command_format = 'bwa mem {ref} {reads_forward} {reads_backward} -o {out}' 
        logger.execute(command_format.format(ref=ref_genome_file,
                                             reads_forward=reads_files[0],
                                             reads_backward=reads_files[1],
                                             out=output_file))
    logger.write('Done.')

    
def compress_sam_file(filename, logger):
    logger.write('-----Compressing {}...'.format(filename))
    # good to check that samtools is intalled
    file_no_ext = os.path.splitext(filename)[0]
    new_file = file_no_ext + '.bam'
    if os.path.exists(new_file):
        logger.write('File {} already exists. Passing'.format(new_file))
    else:
        command_format = 'samtools view -S -b {file_no_ext}.sam > {file_no_ext}.bam'
        logger.execute(command_format.format(file_no_ext=os.path.splitext(filename)[0]))
    logger.write('Done.')


def get_basic_info_aligning(filename, logger):
    logger.write('-----Gathering basic statistics in {}...'.format(filename))
    logger.execute('samtools flagstat {}'.format(filename))
    logger.write('Done.')

    
def sort_and_index_bam_file(filename, logger):
    logger.write('-----Sorting and indexing {}...'.format(filename))
    file_no_ext = os.path.splitext(filename)[0]
    sorted_file = file_no_ext + '_sorted' + '.bam'
    logger.execute('samtools sort {} -o {}'.format(filename,
                                                   sorted_file))
    logger.execute('samtools index {}'.format(sorted_file))
    logger.write('Done.')


def variant_calling(ref_genome_file, sorted_alignment, output_file_pref,
                    varscan_path, logger):
    logger.write('-----Variant calling...')

    command_output, _ = logger.execute('which varscan', can_crash=True)
    if not command_output and not varscan_path:
        logger.write('\n'.join(['varscan not found.',
                                'You need to install it or specify the path to a jar file.',
                                'Cannot proceed. Aborting']))
        exit(0)

    if varscan_path:
        varscan_runner = 'java -jar {}'.format(varscan_path)
    else:
        varscan_runner = 'varscan'
        
    command_format = 'samtools mpileup -f {ref} {alignment} > {out}'
    interm_file = 'my.mpileup'
    if os.path.exists(interm_file):
        logger.write('File {} already exists. Passing'.format(interm_file))
    else:
        logger.execute(command_format.format(ref=ref_genome_file,
                                             alignment=sorted_alignment,
                                             out=interm_file))

    output_file_snp = output_file_pref + '_snp.vcf'
    command_format = ' '.join(['{varscan_runner} mpileup2snp {inp}',
                               '--min-var-freq {min_freq}',
                               '--variants --output-vcf 1',
                               '> {out}'])

    if os.path.exists(output_file_snp):
        logger.write('File {} already exists. Passing'.format(output_file_snp))
    else:
        logger.execute(command_format.format(varscan_runner=varscan_runner,
                                             inp=interm_file,
                                             min_freq=0.2, # need to choose this parameter more carefully
                                             out=output_file_snp))

    output_file_indel = output_file_pref + '_indel.vcf'
    command_format = ' '.join(['{varscan_runner} mpileup2indel {inp}',
                               '--min-var-freq {min_freq}',
                               '--variants --output-vcf 1',
                               '> {out}'])

    if os.path.exists(output_file_indel):
        logger.write('File {} already exists. Passing'.format(output_file_indel))
    else:
        logger.execute(command_format.format(varscan_runner=varscan_runner,
                                             inp=interm_file,
                                             min_freq=0.2, # need to choose this parameter more carefully
                                             out=output_file_indel))

    logger.write('Done')


def get_command_line_args():
    parser = argparse.ArgumentParser(description="Script for finding mutations in DNA")
    parser.add_argument('--trimmomatic', default=None,
                        help='path to the trimmomatic jar file. If not specified, "trimmomatic" command is used')
    parser.add_argument('--varscan', default=None,
                        help='path to the varscan jar file. If not specified, "varscan" command is used')
    args =  parser.parse_args()
    return args


def main():
    args = get_command_line_args()
    raw_data_dir = './raw_data'
    
    if not os.path.exists('files_to_download.csv'):
        print('Cannot find file "files_to_download.csv". Cannot Proceed. Aborting')
        
    # get filenames
    fastq_files = []
    ref_genome_file = None
    with open('files_to_download.csv') as f:
        lines_words = iter(line.strip().split(',') for line in f)
        next(lines_words)
        for _, filename, _ in lines_words:
            ext = os.path.splitext(filename)[1] 
            if ext == '.fastq':
                fastq_files.append(filename)
            if ext == '.fna':
                ref_genome_file = filename
    ref_genome_file_full = os.path.join('./raw_data', ref_genome_file)
        
    logger = Logger()
    get_raw_data(fastq_files + [ref_genome_file],
                 raw_data_dir, logger)

    inspect_raw_data(raw_data_dir, fastq_files, logger)
    inspect_data_fastq(raw_data_dir, fastq_files, logger,
                       '----Inspecting raw data with fastqc...')

    trim_reads(raw_data_dir, fastq_files, logger, args.trimmomatic)
    check_trimmed_reads_count('./trimmomatic_output', logger)
    trimmed_files = ['output_forward_paired.fq', 'output_reverse_paired.fq']
    inspect_data_fastq('./trimmomatic_output', trimmed_files, logger,
                       '-----Inspecting trimmed data with fastqc...')
    index_reference_file(ref_genome_file_full, 
                         logger)
    align_reads(ref_genome_file_full,
                [os.path.join('./trimmomatic_output', filename)
                 for filename in trimmed_files],
                'alignment.sam', logger)
    compress_sam_file('alignment.sam', logger)
    get_basic_info_aligning('alignment.bam', logger)
    sort_and_index_bam_file('alignment.bam', logger)
    variant_calling(ref_genome_file_full,
                    'alignment_sorted.bam',
                    'VarScan_results',
                    args.varscan,
                    logger)
    
    
if __name__ == "__main__":
    main()
