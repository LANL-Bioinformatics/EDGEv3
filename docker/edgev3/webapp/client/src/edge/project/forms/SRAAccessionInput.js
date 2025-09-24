import React, { useState, useEffect } from 'react'
import { isValidSRAAccessionInput } from '../../common/util'
import { TextInput } from './TextInput'
import { components } from './defaults'

export const SRAAccessionInput = (props) => {
  const componentName = 'sraAccessionInput'
  const [form] = useState({ ...components[componentName].init })
  const [validInputs] = useState({ ...components[componentName].validInputs })
  const [doValidation, setDoValidation] = useState(0)

  const setTextInput = (inForm, name) => {
    if (inForm.validForm) {
      form[name] = inForm.textInput.split(/\s*(?:,|$)\s*/)
      form[`${name}_display`] = inForm.textInput
      if (validInputs[name]) {
        validInputs[name].isValid = true
      }
    } else {
      form[name] = []
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
        name={'accessions'}
        setParams={setTextInput}
        defaultValue={props.defaultValue ? props.defaultValue : ''}
        text={props.text ? props.text : components[componentName].params['accessions'].text}
        tooltip={
          props.tooltip ? props.tooltip : components[componentName].params['accessions'].tooltip
        }
        showErrorTooltip={
          props.showErrorTooltip
            ? props.showErrorTooltip
            : components[componentName].params['accessions'].showErrorTooltip
        }
        isOptional={
          props.isOptional
            ? props.isOptional
            : components[componentName].params['accessions'].isOptional
        }
        note={props.note ? props.note : components[componentName].params['accessions'].note}
        placeholder={
          props.placeholder
            ? props.placeholder
            : components[componentName].params['accessions'].placeholder
        }
        errMessage={
          props.errMessage
            ? props.errMessage
            : components[componentName].params['accessions'].errMessage
        }
        isValidTextInput={props.isValidInput ? props.isValidInput : isValidSRAAccessionInput}
      />
    </>
  )
}
