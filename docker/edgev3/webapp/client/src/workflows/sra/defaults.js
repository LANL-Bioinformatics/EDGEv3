import { workflowList } from 'src/util'

export const workflowOptions = [{ value: 'sra2fastq', label: workflowList['sra2fastq'].label }]

export const workflows = {
  sra2fastq: {
    validForm: false,
    errMessage: 'SRA Accession(s) input error',
    inputs: {
      accessions: {
        text: 'SRA Accession(s)',
        value: [],
        display: '',
        sraAccessionInput: {
          tooltip: 'Input SRA accessions (comma separate for > 1 input)',
          placeholder: 'ex: SRR1553609',
          showError: false,
          isOptional: false,
          toUpperCase: true,
          errMessage: 'Invalid SRA accession(s) input',
        },
      },
    },
    validInputs: {
      accessions: { isValid: false, error: 'SRA Accession(s) input error.' },
    },
  },
}
