import React, { useState, useEffect } from 'react'
import { Col, Row } from 'reactstrap'
import MySelect from '../../common/MySelect'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const SelectInput = (props) => {
  const componentName = 'selectInput'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const setNewState2 = (name, value) => {
    setState({
      ...form,
      [name]: value,
    })
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form.selectInput = props.setDefault ? props.options[0].value : ''
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (props.isOptional) {
      form.validForm = true
    } else {
      if (!form.selectInput) {
        form.validForm = false
      }
    }
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`rangeInputTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
            />
          ) : (
            <>{props.text}</>
          )}
        </Col>
        <Col xs="12" md="9">
          <MySelect
            defaultValue={props.setDefault ? props.options[0] : null}
            options={props.options}
            onChange={(e) => {
              if (e) {
                setNewState2('selectInput', e.value)
              } else {
                setNewState2('selectInput', null)
              }
            }}
            placeholder={props.placeholder}
            isClearable={props.isClearable}
          />
        </Col>
      </Row>
    </>
  )
}
