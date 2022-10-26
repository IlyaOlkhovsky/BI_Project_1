import shutil
import subprocess
import os
import argparse


class Logger(object):
    
    def __init__(self, log_file_name):
        self.log_file = log_file_name
        if os.path.exists(log_file_name):
            os.remove(log_file_name)

    def write(self, s, end='\n'):
        print(s, end=end)
        with open(self.log_file, 'a') as f_out:
            f_out.write(s + end)

    def execute(self, command, can_crash=False, target_files=None):
        if target_file and all(map(os.path.exists, target_files)):
            self.write('File {} already exists. Passing')
            return (None, None)
        
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

def check_program_existence(program_name, abort_if_absent=False)
    command_output, _ = logger.execute('which {}'.format(program_name),
                                                         can_crash=True)
    if not command_output and abort_if_absent:
        logger.write(' '.join(['{} not found.'.format(program_name),
                               'You need to install it.',
                               'Cannot proceed. Aborting']))
        exit(0)
        
    return command_output


class Config(object):
    def __init__(self,
                 log_file_name,
                 ref_genome_file,
                 sample_reads_forward_file,
                 sample_reads_backward_file,
                 path_to_fastqc=None,
                 path_to_trimmomatic=None,
                 path_to_varscan=None,
                 path_to_bwa,
                 fastqc_output_dir,
                 trimmomatic_output_dir,
                 bwa_output_dir,
                 trimmomatic_adaptor_file=None):
        
        self.logger = Logger(log_file_name)
        self.ref_genome_file = ref_genome_file
        self.sample_reads_forward_file = sample_reads_forward_file
        self.sample_reads_backward_file = sample_reads_backward_file
        self.path_to_trimmomatic = path_to_trimmomatic
        self.path_to_varscan = path_to_varscan
        self.path_to_bwa = path_to_bwa
        self.fastqc_output_dir = fastqc_output_dir
        self.trimmomatic_output_dir = trimmomatic_output_dir
        self.trimmomatic_adaptor_file = trimmomatic_adaptor_file
        self.bwa_output_dir = bwa_output_dir

        self.check_program_existence()
        self.create_directories()
        self.configure_trimmomatic()
        self.configure_bwa()

    def check_programs_existence(self):
        for prog, path in (('fastqc', self.path_to_fastqc),
                           ('trimmomatic', self.path_to_trimmomatic),
                           ('varscan', self.path_to_varscan),
                           ('bwa', self.path_to_bwa)):
            
            if not path or not os.path.exists(path):
                command_output, _ = self.logger.execute('which {}'.format(program_name),
                                                          can_crash=True)
                if not command_output:
                    self.logger.write('{} not found.'.format(program_name))
                    return False
                
        return True


    def create_directories(self):
        for dir_name in (self.fastqc_output_dir,
                         self.trimmomatic_output_dir,
                         self.bwa_output_dir):
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
                
    def configure_trimmomatic(self):
        trimmed_reads_files = (os.path.splitext(input_file)[0] + '_paired.fq'
                              for input_file in (self.sample_reads_forward_file,
                                                 self.sample_reads_backward_file))
        rest_files = (os.path.splitext(input_file)[0] + '_unpaired.fq'
                      for input_file in (self.sample_reads_forward_file,
                                         self.sample_reads_backward_file))


        self.trimmed_reads_forward_file, self.trimmed_reads_backward_file = trimmed_reads_files
        self.trimmomatic_output_files = trimmed_reads_files + rest_files

        adaptor_path = self.trimmomatic_adaptor_file if self.trimmomatic_adaptor_file \
            else os.path.join(self.path_to_trimmomatic, 'adapters', 'TruSeq3-PE.fa')
        self.trimmomatic_adaptor = 'ILLUMINACLIP:{}:2:30:10'.format(adaptor_path)
    
        if self.path_to_trimmomatic:
            self.trimmomatic_runner = 'java -jar {}'.format(self.path_to_trimmomatic)
        else:
            self.trimmomatic_runner = 'trimmomatic'

           
    def configure_bwa(self):
        shutil.copyfile(self.ref_genome_file,
                        os.path.join(self.bwa_output_dir,
                                     os.path.split(self.ref_genome_file)[1]))
        
def inspect_raw_data(config):
    config.logger.write('----Inspecting raw data...')
    for f in (config.sample_reads_forward_file,
              config.sample_reads_backward_file):
        command_output, _ = config.logger.execute('wc -l {}'.format(f))
        reads_number = int(command_output.split()[0]) // 4
        config.logger.write('Number of reads in {}: {}\n'.format(f, reads_number))

    config.logger.write('Done.\n')

def inspect_data_fastq(files_to_process, config):
    for f in files_to_process:
        command = 'fastqc -o {out} {inp}'.format(out=config.fastqc_output_dir,
                                                 inp=f)
        target_file = os.path.join(config.fastqc_output_dir,
                                   os.path.splitext(f)[0] + '_fastqc.html')
        config.logger.execute(command, target_files=[target_file])

    config.logger.write('Done.\n')
        
def trim_reads(config):
    config.logger.write('-----Trimming reads with Trimmomatic...')

    input_files = ' '.join([config.sample_reads_forward_file,
                            config.sample_reads_backward_file])
    
    output_files = ' '.join(config.trimmomatic_output_files)

    format_string = ' '.join(['{trimmomatic_runner} PE -phred33',
                              '{input_files}',
                              '{output_files}',
                              '{adaptor}',
                              'LEADING:{leading}',
                              'TRAILING:{trailing}',
                              'SLIDINGWINDOW:{sliding_window_1}:{sliding_window_2}',
                              'MINLEN:{minlen}'])

    config.logger.execute(format_string.format(trimmomatic_runner=config.trimmomatic_runner,
                                               input_files=input_files,
                                               output_files=output_files,
                                               adaptor=config.trimmomatic_adaptor,
                                               leading=20,
                                               trailing=20,
                                               sliding_window_1=10,
                                               sliding_window_2=20,
                                               minlen=20),
                          target_files=output_files.split())
    
    config.logger.write('Done.\n')
                                                                                                                                                                                                    
def check_trimmed_reads_count(config):
    config.logger.write('-----Checking trimmed reads count...')
    for filename in (config.trimmed_reads_forward_file,
                     config.trimmed_reads_backward_file):
        config.logger.execute('zcat {} | wc -l'.format(filename))

def index_reference_file(config):
    logger.write('-----Indexing reference file...')
    
    ref_file = os.path.join(config.bwa_output_dir,
                            os.path.split(config.ref_genome_file)[1])
    
    shutil.copyfile(config.ref_genome_file, ref_file)
    
    logger.execute('bwa index {}'.format(ref_file))
    logger.write('Done.\n')



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
    logger.write('Done.\n')

    
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
    download_list = []
    fastq_files = []
    ref_genome_file = None
    with open('files_to_download.csv') as f:
        lines_words = iter(line.strip().split(',') for line in f)
        next(lines_words)
        for download_link, filename, _ in lines_words:
            download_list.append((filename, download_link))
            ext = os.path.splitext(filename)[1] 
            if ext == '.fastq':
                fastq_files.append(filename)
            if ext == '.fna':
                ref_genome_file = filename
    ref_genome_file_full = os.path.join('./raw_data', ref_genome_file)
        
    logger = Logger()
    get_raw_data(download_list,
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
