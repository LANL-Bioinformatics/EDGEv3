import React, { useState, useEffect } from 'react'
import { Col, Row } from 'reactstrap'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'
import { RangeTextInput } from './RangeTextInput'

export const RangeTextInputSingle = (props) => {
  const componentName = 'rangeTextInputSingle'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const setRangeTextInputSingle = (params, name) => {
    setState({ ...params })
    setDoValidation(doValidation + 1)
  }

  //trigger validation method when input changes
  useEffect(() => {
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`rangeTextInputSingleTooltip-${props.name}`}
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
          <RangeTextInput
            name={'rangeTextInputSingle'}
            setParams={setRangeTextInputSingle}
            startText={props.startText}
            endText={props.endText}
            defaultValueStart={props.defaultValueStart}
            minStart={props.minStart}
            maxStart={props.maxStart}
            defaultValueEnd={props.defaultValueEnd}
            minEnd={props.minEnd}
            maxEnd={props.maxEnd}
          />
        </Col>
      </Row>
    </>
  )
}
