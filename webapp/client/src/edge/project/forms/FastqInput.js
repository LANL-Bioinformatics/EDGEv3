import React, { useState, useEffect } from 'react'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { Switcher } from './Switcher'
import { FileInputArray } from './FileInputArray'
import { PairedFileInputArray } from './PairedFileInputArray'
import { components } from './defaults'

export const FastqInput = (props) => {
  const componentName = 'fastqInput'
  const [form, setState] = useState({ ...components[componentName].init })
  const [doValidation, setDoValidation] = useState(0)

  const setNewState2 = (name, value) => {
    setState({
      ...form,
      [name]: value,
    })
    setDoValidation(doValidation + 1)
  }
  const setSwitcher = (inForm, name) => {
    setNewState2(name, inForm.isTrue)
  }
  const setFileInput = (inForm, name) => {
    form.validForm = inForm.validForm
    if (inForm.validForm) {
      setState({
        ...form,
        fileInput: inForm.fileInput,
        fileInput_display: inForm.fileInput_display,
      })
    } else {
      setState({
        ...form,
        fileInput: [],
        fileInput_display: [],
      })
    }
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    //set interleaved
    if (props.isInterleaved != null) {
      setNewState2('interleaved', props.isInterleaved)
    }
  }, [props.isInterleaved]) // eslint-disable-line react-hooks/exhaustive-deps

  //trigger validation method when input changes
  useEffect(() => {
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      {props.text && (
        <MyTooltip
          id={`fileInputTooltip-${props.name}`}
          tooltip={props.tooltip}
          text={props.text}
          place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
          color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
          showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
        />
      )}
      {!props.disableSwitcher && (
        <>
          <Switcher
            id={'interleaved'}
            name={'interleaved'}
            setParams={setSwitcher}
            text={components[componentName].params['interleaved'].text}
            defaultValue={components[componentName].params['interleaved'].defaultValue}
            trueText={components[componentName].params['interleaved'].trueText}
            falseText={components[componentName].params['interleaved'].falseText}
          />
          <br></br>
        </>
      )}

      {form.interleaved && (
        <>
          <FileInputArray
            setParams={setFileInput}
            name={'fastq'}
            text={components[componentName].params['fastq'].text}
            enableInput={props.enableInput}
            placeholder={props.placeholder}
            isValidFileInput={props.isValidFileInput}
            dataSources={props.dataSources}
            fileTypes={props.fileTypes}
            projectTypes={props.projectTypes}
            viewFile={props.viewFile}
            isOptional={props.isOptional}
            cleanupInput={props.cleanupInput}
            maxInput={props.maxInput}
          />
        </>
      )}
      {!form.interleaved && (
        <>
          <PairedFileInputArray
            setParams={setFileInput}
            name={'fastq'}
            text={components[componentName].params['pairedFastq'].text}
            enableInput={props.enableInput}
            placeholder={props.placeholder}
            isValidFileInput={props.isValidFileInput}
            dataSources={props.dataSources}
            fileTypes={props.fileTypes}
            projectTypes={props.projectTypes}
            viewFile={props.viewFile}
            isOptional={props.isOptional}
            cleanupInput={props.cleanupInput}
            maxInput={props.maxInput}
          />
        </>
      )}
    </>
  )
}
