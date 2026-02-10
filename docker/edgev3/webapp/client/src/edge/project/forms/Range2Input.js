import React, { useState, useEffect } from 'react'
import { Col, Row } from 'reactstrap'
import Slider, { Range } from 'rc-slider'
import 'rc-slider/assets/index.css'
import 'rc-slider/assets/index.css'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const Range2Input = (props) => {
  const componentName = 'range2Input'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const setRange2 = (value) => {
    form['start'] = value[0]
    form['end'] = value[1]
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form['start'] = props.defaultValue[0]
    form['end'] = props.defaultValue[1]
  }, [props.defaultValue]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`range2InputTooltip-${props.name}`}
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
          <Row>
            <Col> {form['start']}</Col>
            <Col className="edge-right"> {form['end']} </Col>
          </Row>

          <Slider
            range
            name="range2Input"
            value={[form['start'], form['end']]}
            min={props.min}
            max={props.max}
            step={props.step}
            onChange={(e) => setRange2(e)}
          />
        </Col>
      </Row>
    </>
  )
}
