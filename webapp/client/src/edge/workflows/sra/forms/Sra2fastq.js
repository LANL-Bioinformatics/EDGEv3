import React, { useState, useEffect } from 'react'
import { Card, CardBody, Collapse } from 'reactstrap'
import { Header } from '../../../project/forms/SectionHeader'
import { TextInput } from '../../../project/forms/TextInput'
import { workflows } from '../defaults'

export const Sra2fastq = (props) => {
  const workflowName = 'sra2fastq'
  const [collapseParms, setCollapseParms] = useState(false)
  const [form] = useState({ ...workflows[workflowName] })
  const [validInputs] = useState({ ...workflows[workflowName].validInputs })
  const [doValidation, setDoValidation] = useState(0)

  const toggleParms = () => {
    setCollapseParms(!collapseParms)
  }

  const isValidAccessions = (accessions) => {
    if (!accessions) {
      if (props.required) {
        return false
      } else {
        return true
      }
    }
    const parts = accessions.split(',')
    for (var i = 0; i < parts.length; i++) {
      //if(!/^[a-zA-Z]{3}[0-9]{6,9}$/.test(parts[i].trim())) {
      if (
        !/^(srp|erp|drp|srx|erx|drx|srs|ers|drs|srr|err|drr|sra|era|dra)[0-9]{6,9}$/i.test(
          parts[i].trim(),
        )
      ) {
        return false
      }
    }
    return true
  }

  const setTextInput = (inForm, name) => {
    if (inForm.validForm) {
      form.inputs[name].value = inForm.textInput.split(/\s*(?:,|$)\s*/)
      form.inputs[name].display = inForm.textInput
      if (validInputs[name]) {
        validInputs[name].isValid = true
      }
    } else {
      form.inputs[name].value = []
      if (validInputs[name]) {
        validInputs[name].isValid = false
      }
    }
    setDoValidation(doValidation + 1)
  }

  //trigger validation method when input changes
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
    <Card className="workflow-card">
      <Header
        toggle={true}
        toggleParms={toggleParms}
        title={'Input'}
        collapseParms={collapseParms}
        id={workflowName + 'input'}
        isValid={props.isValid}
        errMessage={props.errMessage}
      />
      <Collapse isOpen={!collapseParms} id={'collapseParameters-' + props.name}>
        <CardBody>
          <TextInput
            name={'accessions'}
            setParams={setTextInput}
            text={workflows[workflowName].inputs['accessions'].text}
            tooltip={workflows[workflowName].inputs['accessions']['textInput'].tooltip}
            showError={workflows[workflowName].inputs['accessions']['textInput'].showError}
            isOptional={workflows[workflowName].inputs['accessions']['textInput'].isOptional}
            toUpperCase={workflows[workflowName].inputs['accessions']['textInput'].toUpperCase}
            note={workflows[workflowName].inputs['accessions']['textInput'].note}
            placeholder={workflows[workflowName].inputs['accessions']['textInput'].placeholder}
            errMessage={workflows[workflowName].inputs['accessions']['textInput'].errMessage}
            isValidTextInput={isValidAccessions}
            reset={props.reset}
          />
        </CardBody>
      </Collapse>
    </Card>
  )
}
