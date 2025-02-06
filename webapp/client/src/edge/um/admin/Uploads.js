import React from 'react'
import { Badge } from 'reactstrap'
import UploadTable from '../common/UploadTable'

const Uploads = (props) => {
  return (
    <div className="animated fadeIn">
      <div className="clearfix">
        <Badge color="danger" pill>
          Admin tool
        </Badge>
      </div>
      <br></br>
      <UploadTable tableType="admin" title={'Manage Uploads'} {...props} />
      <br></br>
    </div>
  )
}

export default Uploads
