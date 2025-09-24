import React, { useState } from 'react'
import { Button, ButtonGroup, CardHeader } from 'reactstrap'
import CIcon from '@coreui/icons-react'
import { cilChevronTop, cilChevronBottom } from '@coreui/icons'
import { colors } from '../../common/util'

export const Header = (props) => {
  const [headerColor, setHeaderColor] = useState(colors.light)

  return (
    <CardHeader
      style={{ backgroundColor: headerColor, cursor: 'pointer' }}
      onMouseEnter={(e) => setHeaderColor(colors.secondary)}
      onMouseLeave={() => setHeaderColor(colors.light)}
      onClick={props.toggleParms}
    >
      {props.toggle ? (
        <>
          <span
            className="edge-card-header-action"
            data-target="#collapseParameters"
            onClick={props.toggleParms}
          >
            {!props.collapseParms ? (
              <CIcon icon={cilChevronTop} />
            ) : (
              <CIcon icon={cilChevronBottom} />
            )}
          </span>
          &nbsp;&nbsp;&nbsp;
        </>
      ) : (
        <>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</>
      )}
      <span className="edge-card-header">{props.title}</span>
      {props.onoff && (
        <ButtonGroup style={{ float: 'right' }} className="mr-3" aria-label="First group" size="sm">
          <Button
            color="outline-primary"
            onClick={() => props.setOnoff(true)}
            active={props.paramsOn}
          >
            On
          </Button>
          <Button
            color="outline-primary"
            onClick={() => {
              props.setOnoff(false)
              props.setCollapseParms(true)
            }}
            active={!props.paramsOn}
          >
            Off
          </Button>
        </ButtonGroup>
      )}
    </CardHeader>
  )
}
