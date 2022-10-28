import sys
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
        if target_files and all(map(os.path.exists, target_files)):
            self.write('Files {} already exist. Passing'.format(target_files))
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


def file_name_no_ext(full_name):
    return os.path.splitext(os.path.split(full_name)[1])[0]


class Config(object):

    def check_programs_existence(self):
        for prog, path in (('fastqc', self.path_to_fastqc),
                           ('trimmomatic', self.path_to_trimmomatic),
                           ('varscan', self.path_to_varscan),
                           ('bwa', self.path_to_bwa),
                           ('samtools', None)):

            if not path or not os.path.exists(path):
                command_output, _ = self.logger.execute('which {}'.format(prog),
                                                        can_crash=True)
                if not command_output:
                    self.logger.write('{} not found.'.format(prog))
                    return False

        return True

    def create_directories(self):
        for dir_name in (self.fastqc_output_dir,
                         self.trimmomatic_output_dir,
                         self.bwa_output_dir,
                         self.varscan_output_dir):
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)

    def configure_trimmomatic(self):
        self.trimmed_reads_forward_file, self.trimmed_reads_backward_file = \
            [os.path.join(self.trimmomatic_output_dir,
                          file_name_no_ext(input_file) + '_paired.fq')
             for input_file in (self.sample_reads_forward_file,
                                self.sample_reads_backward_file)]

        self.trimmed_reads_unpaired_forward, self.trimmed_reads_unpaired_backward = \
            [os.path.join(self.trimmomatic_output_dir,
                          file_name_no_ext(input_file) + '_unpaired.fq')
             for input_file in (self.sample_reads_forward_file,
                                self.sample_reads_backward_file)]

        adaptor_path = self.trimmomatic_adaptor_file if self.trimmomatic_adaptor_file \
            else os.path.join(os.path.split(self.path_to_trimmomatic)[0], 'adapters', 'TruSeq3-PE.fa')
        self.trimmomatic_adaptor = 'ILLUMINACLIP:{}:2:30:10'.format(adaptor_path)

        if self.path_to_trimmomatic:
            self.trimmomatic_runner = 'java -jar {}'.format(self.path_to_trimmomatic)
        else:
            self.trimmomatic_runner = 'trimmomatic'


    def configure_bwa(self):
        self.alignment_file_compressed = \
            os.path.splitext(self.alignment_file)[0] + '.bam'
        self.alignment_file_sorted = \
            os.path.splitext(self.alignment_file_compressed)[0] + \
            '_sorted' + '.bam'

    def configure_varscan(self):
        if self.path_to_varscan:
            self.varscan_runner = 'java -jar {}'.format(self.path_to_varscan)
        else:
            self.varscan_runner = 'varscan'


    def __init__(self, *,
                 log_file_name,
                 ref_genome_file,
                 sample_reads_forward_file,
                 sample_reads_backward_file,
                 path_to_fastqc=None,
                 path_to_trimmomatic=None,
                 path_to_varscan=None,
                 path_to_bwa=None,
                 fastqc_output_dir,
                 trimmomatic_output_dir,
                 bwa_output_dir,
                 varscan_output_dir,
                 trimmomatic_adaptor_file=None,
                 alignment_file,
                 varscan_interm_file,
                 varscan_result_pref):

        self.logger = Logger(log_file_name)
        self.ref_genome_file = ref_genome_file
        self.sample_reads_forward_file = sample_reads_forward_file
        self.sample_reads_backward_file = sample_reads_backward_file
        self.path_to_fastqc = path_to_fastqc
        self.path_to_trimmomatic = path_to_trimmomatic
        self.path_to_varscan = path_to_varscan
        self.path_to_bwa = path_to_bwa
        self.fastqc_output_dir = fastqc_output_dir
        self.trimmomatic_output_dir = trimmomatic_output_dir
        self.varscan_output_dir = varscan_output_dir
        self.trimmomatic_adaptor_file = trimmomatic_adaptor_file
        self.bwa_output_dir = bwa_output_dir
        self.alignment_file = os.path.join(bwa_output_dir,
                                           alignment_file)
        self.varscan_interm_file = os.path.join(varscan_output_dir,
                                                varscan_interm_file)
        self.varscan_result_pref = os.path.join(varscan_output_dir,
                                                varscan_result_pref)

        if not self.check_programs_existence():
            self.logger.write('Cannot proceed. Aborting')
            exit(0)

        self.create_directories()
        self.configure_trimmomatic()
        self.configure_bwa()
        self.configure_varscan()


def inspect_raw_data(config):
    for f in (config.sample_reads_forward_file,
              config.sample_reads_backward_file):
        command_output, _ = config.logger.execute('wc -l {}'.format(f))
        reads_number = int(command_output.split()[0]) // 4
        config.logger.write('Number of reads in {}: {}\n'.format(f, reads_number))


def inspect_data_fastq(files_to_process, config):
    command = 'fastqc -o {out} {inp}'.format(out=config.fastqc_output_dir,
                                             inp=' '.join(files_to_process))
    target_files = [os.path.join(config.fastqc_output_dir,
                                 file_name_no_ext(f) + '_fastqc.html')
                    for f in files_to_process]

    config.logger.execute(command, target_files=target_files)


def trim_reads(config):
    input_files = ' '.join([config.sample_reads_forward_file,
                            config.sample_reads_backward_file])

    output_files = ' '.join([config.trimmed_reads_forward_file,
                             config.trimmed_reads_unpaired_forward,
                             config.trimmed_reads_backward_file,
                             config.trimmed_reads_unpaired_backward])

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


def check_trimmed_reads_count(config):
    for filename in (config.trimmed_reads_forward_file,
                     config.trimmed_reads_backward_file):
        config.logger.execute('cat {} | wc -l'.format(filename))


def index_reference_file(config):
    ref_file = os.path.join(config.bwa_output_dir,
                            os.path.split(config.ref_genome_file)[1])

    shutil.copyfile(config.ref_genome_file, ref_file)
    config.logger.execute('bwa index {}'.format(ref_file),
                          target_files=[ref_file + ext
                                        for ext in ('.amb', '.ann', '.bwt', '.pac', '.sa')])


def align_reads(config):
    command_format = 'bwa mem {ref} {reads_forward} {reads_backward} -o {out}'
    ref_file = os.path.join(config.bwa_output_dir,
                            os.path.split(config.ref_genome_file)[1])
    config.logger.execute(command_format.format(ref=ref_file,
                                                reads_forward=config.trimmed_reads_forward_file,
                                                reads_backward=config.trimmed_reads_backward_file,
                                                out=config.alignment_file),
                          target_files=[config.alignment_file])


def compress_sam_file(config):
    command_format = 'samtools view -S -b {inp} > {out}'
    config.logger.execute(command_format.format(inp=config.alignment_file,
                                                out=config.alignment_file_compressed),
                          target_files=[config.alignment_file_compressed])


def get_basic_info_aligning(config):
    config.logger.execute('samtools flagstat {}'.format(config.alignment_file_compressed))


def sort_and_index_bam_file(config):
    config.logger.execute('samtools sort {} -o {}'.format(config.alignment_file_compressed,
                                                          config.alignment_file_sorted),
                          target_files=[config.alignment_file_sorted])
    config.logger.execute('samtools index {}'.format(config.alignment_file_sorted),
                          target_files=[config.alignment_file_sorted + '.bai'])


def variant_calling(config):
    command_format = 'samtools mpileup -f {ref} {alignment} > {out}'
    config.logger.execute(command_format.format(ref=config.ref_genome_file,
                                                alignment=config.alignment_file_sorted,
                                                out=config.varscan_interm_file),
                          target_files=[config.varscan_interm_file])

    output_file_snp = config.varscan_result_pref + '_snp.vcf'
    command_format = ' '.join(['{varscan_runner} mpileup2snp {inp}',
                               '--min-var-freq {min_freq}',
                               '--variants --output-vcf 1',
                               '> {out}'])

    config.logger.execute(command_format.format(varscan_runner=config.varscan_runner,
                                                inp=config.varscan_interm_file,
                                                min_freq=0.2, # need to choose this parameter more carefully
                                                out=output_file_snp),
                          target_files=[output_file_snp])

    output_file_indel = config.varscan_result_pref + '_indel.vcf'
    command_format = ' '.join(['{varscan_runner} mpileup2indel {inp}',
                               '--min-var-freq {min_freq}',
                               '--variants --output-vcf 1',
                               '> {out}'])

    config.logger.execute(command_format.format(varscan_runner=config.varscan_runner,
                                                inp=config.varscan_interm_file,
                                                min_freq=0.2, # need to choose this parameter more carefully
                                                out=output_file_indel),
                          target_files=[output_file_indel])


def get_command_line_args():
    parser = argparse.ArgumentParser(description="Script for finding mutations in DNA")
    parser.add_argument('--trimmomatic', default=None,
                        help='path to the trimmomatic jar file. If not specified, "trimmomatic" command is used')
    parser.add_argument('--varscan', default=None,
                        help='path to the varscan jar file. If not specified, "varscan" command is used')
    parser.add_argument('--ref_genome',
                        help='path to the .fna file with reference genome')
    parser.add_argument('--reads_forward',
                        help='path to the .fastq file with forward read reads')
    parser.add_argument('--reads_reverse',
                        help='path to the .fastq file with reverse read reads')

    args = parser.parse_args()
    return args


def main():
    print('Started BI_Project_1 with the arguments\n')
    print('"' + ' '.join(sys.argv[1:]) + '"')

    args = get_command_line_args()
    config = Config(log_file_name='log.txt',
                    ref_genome_file=args.ref_genome,
                    sample_reads_forward_file=args.reads_forward,
                    sample_reads_backward_file=args.reads_reverse,
                    path_to_trimmomatic=args.trimmomatic,
                    path_to_varscan=args.varscan,
                    fastqc_output_dir='./fastqc_output',
                    trimmomatic_output_dir='./trimmomatic_output',
                    bwa_output_dir='./bwa_output',
                    alignment_file='alignment.sam',
                    varscan_output_dir='./varscan_output',
                    varscan_interm_file='my.mpileup',
                    varscan_result_pref='VarScan_results')

    config.logger.write('-----Inspecting raw data...')
    inspect_raw_data(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Processing raw data with fastqc...')
    inspect_data_fastq(files_to_process=[config.sample_reads_forward_file,
                                         config.sample_reads_backward_file],
                       config=config)
    config.logger.write('Done.\n')

    config.logger.write('-----Trimming data to get rid of low-quality reads...')
    trim_reads(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Checking trimmed reads count...')
    check_trimmed_reads_count(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Processing trimmed data with fastqc...')
    inspect_data_fastq(files_to_process=[config.trimmed_reads_forward_file,
                                         config.trimmed_reads_backward_file],
                       config=config)
    config.logger.write('Done.\n')

    config.logger.write('-----Indexing reference file...')
    index_reference_file(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Aligning reads...')
    align_reads(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Compressing alignment file...')
    compress_sam_file(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Gathering basic statistics in alignment file...')
    get_basic_info_aligning(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Sorting and indexing compressed alignment file...')
    sort_and_index_bam_file(config)
    config.logger.write('Done.\n')

    config.logger.write('-----Variant calling...')
    variant_calling(config)
    config.logger.write('Done.\n')


if __name__ == "__main__":
    main()
