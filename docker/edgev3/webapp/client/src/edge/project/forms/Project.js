import React, { useState, useEffect } from 'react'
import { isValidProjectName } from '../../common/util'
import { TextInput } from './TextInput'
import { components } from './defaults'

export const Project = (props) => {
  const componentName = 'project'
  const [form, setState] = useState({ ...components[componentName].init })
  const [validInputs] = useState({ ...components[componentName].validInputs })
  const [doValidation, setDoValidation] = useState(0)

  const setTextInput = (inForm, name) => {
    if (inForm.validForm) {
      setState({
        ...form,
        [name]: inForm.textInput,
      })
      if (validInputs[name]) {
        validInputs[name].isValid = true
      }
    } else {
      setState({
        ...form,
        [name]: null,
      })
      if (validInputs[name]) {
        validInputs[name].isValid = false
      }
    }
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    // check input errors
    let errors = ''
    Object.keys(validInputs).forEach((key) => {
      if (!validInputs[key].isValid) {
        errors += validInputs[key].error + '<br/>'
      }
    })

    if (errors === '') {
      form.errMessage = null
      form.validForm = true
    } else {
      form.errMessage = errors
      form.validForm = false
    }
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <TextInput
        name={'projectName'}
        setParams={setTextInput}
        defaultValue={''}
        text={components[componentName].params['projectName'].text}
        showErrorTooltip={components[componentName].params['projectName'].showErrorTooltip}
        isOptional={components[componentName].params['projectName'].isOptional}
        note={components[componentName].params['projectName'].note}
        placeholder={components[componentName].params['projectName'].placeholder}
        errMessage={components[componentName].params['projectName'].errMessage}
        isValidTextInput={isValidProjectName}
      />
      <br></br>
      <TextInput
        name={'projectDesc'}
        setParams={setTextInput}
        defaultValue={''}
        text={components[componentName].params['projectDesc'].text}
        isOptional={components[componentName].params['projectDesc'].isOptional}
        placeholder={components[componentName].params['projectDesc'].placeholder}
        isValidTextInput={() => {
          return true
        }}
      />
    </>
  )
}
