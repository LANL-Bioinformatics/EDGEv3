import React, { useState } from 'react'
import UploadTable from '../common/UploadTable'
import UploadFiles from './UploadFiles'

const Uploads = (props) => {
  const [refreshTable, setRefreshTable] = useState(0)

  const reloadTableData = () => {
    setRefreshTable(refreshTable + 1)
  }
  return (
    <div className="animated fadeIn">
      <UploadFiles reloadTableData={reloadTableData} />
      <br></br>
      <br></br>
      <UploadTable tableType={'user'} title={'My Uploads'} refresh={refreshTable} />
      <br></br>
    </div>
  )
}

export default Uploads
