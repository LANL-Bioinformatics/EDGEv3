import React from 'react'

const HtmlToReactParser = require('html-to-react').Parser
let htmlToReactParser = new HtmlToReactParser()

export const HtmlText = (props) => {
  return <>{props.text && <>{htmlToReactParser.parse(props.text)}</>}</>
}
