export const components = {
  textInput: {
    validForm: false,
    textInput: '',
  },
  textInputArray: {
    validForm: false,
    textInputs: [],
  },
  integerInput: {
    validForm: false,
    integerInput: null,
  },
  rangeInput: {
    validForm: true,
    rangeInput: 0,
  },
  rangeTextInput: {
    validForm: false,
    start: 0,
    end: 1,
  },
  rangeTextInputArray: {
    validForm: false,
    rangeInputs: [],
  },
  switcher: {
    validForm: true,
    isTrue: true,
  },
  selectInput: {
    validForm: false,
    selectInput: null,
  },
  fileUpload: {
    file: '',
    validForm: false,
    errMessage: 'File is required',
  },
  fileInput: {
    validForm: false,
    fileInput: null,
    fileInput_display: null,
  },
  fileInputArray: {
    validForm: false,
    fileInput: [],
    fileInput_display: [],
    fileInput_isValid: [],
  },
  pairedFileInputArray: {
    validForm: false,
    fileInput: [],
    fileInput_display: [],
    fileInput_isValid: [],
  },
  fastqInput: {
    params: {
      interleaved: {
        trueText: 'Yes',
        falseText: 'No',
        defaultValue: true,
        text: 'Is interleaved?',
      },
      fastq: {
        text: 'Fastq',
      },
      pairedFastq: {
        text: 'Paired Fastq',
      },
    },
    init: {
      validForm: false,
      errMessage: '',
      interleaved: true,
      fileInput: [],
      fileInput_display: [],
    },
  },
  project: {
    params: {
      projectName: {
        text: 'Project/Run Name',
        placeholder: 'required, at 3 but less than 30 characters',
        showError: false,
        isOptional: false,
        showErrorTooltip: true,
        errMessage:
          'Required, at 3 but less than 30 characters. <br/>Only alphabets, numbers, dashs, dot and underscore are allowed in the name.',
      },
      projectDesc: {
        text: 'Description',
        placeholder: 'optional',
        showError: false,
        isOptional: true,
        showErrorTooltip: false,
        errMessage: '',
      },
    },
    validInputs: {
      projectName: { isValid: false, error: 'Project/Run Name input error.' },
    },
    init: {
      validForm: false,
      errMessage: null,
      projectName: null,
      projectDesc: null,
    },
  },
}
