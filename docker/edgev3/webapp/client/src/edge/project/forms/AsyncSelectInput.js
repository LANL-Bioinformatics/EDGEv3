import React, { useState, useEffect } from 'react'
import { Col, Row } from 'reactstrap'
import { MyAsyncSelect } from '../../common/MySelect'
import { MyTooltip, ErrorTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const AsyncSelectInput = (props) => {
  const componentName = 'asyncSelectInput'
  const [form, setState] = useState({ ...components[componentName] })
  const [noOptionMsg, setNoOptionMsg] = useState('No result')

  const [doValidation, setDoValidation] = useState(0)

  useEffect(() => {
    form.selections = props.value ? props.value : []
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  const handleSelectChange = (selected) => {
    if (selected) {
      form.selections = selected
    } else {
      form.selections = []
    }
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form.selections = props.value ? props.value : []
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    form.validForm = true
    form.errMessage = null
    const min = props.min ? props.min : 1
    const max = props.max ? props.max : 1000000
    if (!props.isOptional && (form.selections.length < min || form.selections.length > max)) {
      form.validForm = false
      form.errMessage = `Select at least ${min} but no more than ${max} items`
    }
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`multSelctTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
            />
          ) : (
            <>
              {props.text}
              {errors && errors.textInput && props.showErrorTooltip && (
                <ErrorTooltip
                  id={`textInputErrTooltip-${props.name}`}
                  tooltip={errors.textInput.message}
                />
              )}
            </>
          )}
          <br></br>
          {props.showSelections && form.selections.length > 0 && (
            <>
              <span className="text-muted edge-text-size-small">
                {form.selections.length} {props.showSelectionsText}
              </span>
              <br></br>
            </>
          )}
        </Col>
        <Col xs="12" md="9">
          <MyAsyncSelect
            options={props.options}
            value={props.value}
            isMulti={props.isMulti}
            placeholder={props.placeholder ? props.placeholder : 'Enter term(s) to search...'}
            closeMenuOnSelect={false}
            hideSelectedOptions={false}
            onChange={handleSelectChange}
            isClearable={true}
            selectAll={true}
            maxOptions={props.maxOptions}
            maxSelected={props.max}
            setMsg={(msg) => setNoOptionMsg(msg)}
            noOptionsMessage={noOptionMsg}
          />
        </Col>
      </Row>
      <br></br>
    </>
  )
}
