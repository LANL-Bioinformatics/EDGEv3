import { workflowList } from 'src/util'

export const workflowOptions = [{ value: 'runFaQCs', label: workflowList['runFaQCs'].label }]

export const workflows = {
  runFaQCs: {
    validForm: false,
    errMessage: 'input error',
    inputs: {
      interleaved: {
        text: 'Is interleaved',
        value: true,
      },
      inputFastq: {
        text: 'Input Raw Reads',
        value: [],
        display: [],
        fastqInput: {
          tooltip:
            'RunFaQCs requires paired-end Illumina data in FASTQ format as the input; the file can be interleaved and can becompressed. <br/>Acceptable file formats: .fastq, .fq, .fastq.gz, .fq.gz',
          enableInput: true,
          placeholder: 'Select a file or enter a file http(s) url',
          dataSources: ['upload', 'public', 'project'],
          fileTypes: ['fastq', 'fq', 'fastq.gz', 'fq.gz'],
          projectTypes: ['sra2fastq'],
          viewFile: false,
          isOptional: false,
          cleanupInput: true,
          maxInput: 1,
        },
      },
      trimQual: {
        text: 'Trim Quality Level',
        value: 20,
        rangeInput: {
          tooltip: 'Targets # as quality level (default 20) for trimming',
          defaultValue: 20,
          min: 0,
          max: 40,
          step: 1,
        },
      },
      trim5end: {
        text: "Cut #bp from 5'-end",
        value: 0,
        rangeInput: {
          tooltip: 'Cut # bp from 5 end before quality trimming/filtering',
          defaultValue: 0,
          min: 0,
          max: 100,
          step: 1,
        },
      },
      trim3end: {
        text: "Cut #bp from 3'-end",
        value: 0,
        rangeInput: {
          tooltip: 'Cut # bp from 3 end before quality trimming/filtering',
          defaultValue: 0,
          min: 0,
          max: 100,
          step: 1,
        },
      },
      trimAdapter: {
        text: ' Trim Adapter',
        value: false,
        switcher: {
          tooltip: 'Trim reads with illumina adapter/primers (default: No)',
          trueText: 'Yes',
          falseText: 'No',
          defaultValue: false,
        },
      },
      trimRate: {
        text: 'Trim Adapter mismatch ratio',
        value: 0.2,
        rangeInput: {
          tooltip: "Mismatch ratio of adapters' length(default : 0.2, allow 20% mismatches) ",
          defaultValue: 0.2,
          min: 0.0,
          max: 1.0,
          step: 0.01,
        },
      },
      trimPolyA: {
        text: 'Trim polyA',
        value: false,
        switcher: {
          tooltip: 'Trim poly A ( > 15 )',
          trueText: 'Yes',
          falseText: 'No',
          defaultValue: false,
        },
      },
      artifactFile: {
        text: 'Adapter/Primer FASTA',
        value: null,
        display: null,
        fileInput: {
          tooltip:
            'Additional artifact (adapters/primers/contaminations) reference file in fasta format',
          enableInput: true,
          placeholder: '(Optional) Select a file or enter a file http(s) url',
          dataSources: ['upload', 'public'],
          fileTypes: ['fasta', 'fa', 'fna', 'contigs'],
          viewFile: false,
          isOptional: true,
          cleanupInput: true,
        },
      },
      minLen: {
        text: 'Minimum Read Length',
        value: 50,
        integerInput: {
          tooltip:
            'Trimmed read should have to be at least this minimum length (default:50, range: 0 - 1000)',
          defaultValue: 50,
          min: 0,
          max: 1000,
        },
      },
      avgQual: {
        text: 'Average Quality Cutoff',
        value: 0,
        rangeInput: {
          tooltip: 'Average quality cutoff (default:0, no filtering)',
          defaultValue: 0,
          min: 0,
          max: 40,
          step: 1,
        },
      },
      numN: {
        text: '"N" Base Cutoff',
        value: 2,
        rangeInput: {
          tooltip:
            'Trimmed read has greater than or equal to this number of continuous base "N" will be discarded. (default: 2, "NN")',
          defaultValue: 2,
          min: 1,
          max: 10,
          step: 1,
        },
      },
      filtLC: {
        text: 'Low Complexity Filter',
        value: 0.85,
        rangeInput: {
          tooltip:
            'Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)',
          defaultValue: 0.85,
          min: 0.0,
          max: 1.0,
          step: 0.01,
        },
      },
      filtPhiX: {
        text: 'Filter phiX',
        value: false,
        switcher: {
          tooltip: 'Filter phiX reads',
          trueText: 'Yes',
          falseText: 'No',
          defaultValue: false,
        },
      },
    },
    // only for input with validation method
    validInputs: {
      inputFastq: { isValid: false, error: 'Input Raw Reads error.' },
      artifactFile: { isValid: true, error: 'Adapter/Primer FASTA error.' },
      minLen: { isValid: true, error: 'Minimum Read Length error.' },
    },
  },
}
