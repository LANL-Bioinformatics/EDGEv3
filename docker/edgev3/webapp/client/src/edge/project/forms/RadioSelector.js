import React, { useState, useEffect } from 'react'
import { Col, Row, Input, FormGroup, Label } from 'reactstrap'
import { HtmlText } from '../../common/HtmlText'

import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const RadioSelector = (props) => {
  const componentName = 'radioSelector'
  const [form, setState] = useState({
    ...components[componentName],
    option: props.defaultValue,
    display: props.display ? props.display : props.defaultValue,
  })
  const [doValidation, setDoValidation] = useState(0)

  useEffect(() => {
    setState({
      ...components[componentName],
      option: props.defaultValue,
      display: props.display ? props.display : props.defaultValue,
    })
  }, [props.defaultValue]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    setState({
      ...components[componentName],
      option: props.defaultValue,
      display: props.display ? props.display : props.defaultValue,
    })
  }, [props.reset]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md={props.text ? '3' : '0'}>
          {props.tooltip ? (
            <MyTooltip
              id={`radioSelectorTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
            />
          ) : (
            <>{props.text}</>
          )}{' '}
        </Col>
        <Col xs="12" md={props.text ? '9' : '12'}>
          {props.options.map((item, index) => (
            <FormGroup check key={`group-${index}`}>
              <Input
                key={`input-${index}`}
                type="radio"
                checked={form.option === item.value}
                onChange={() => {
                  form.display = item.text
                  form.option = item.value
                  setDoValidation(doValidation + 1)
                }}
              />{' '}
              <Label check key={`label-${index}`}>
                <HtmlText text={item.detail} />
              </Label>
            </FormGroup>
          ))}
        </Col>
      </Row>
    </>
  )
}
