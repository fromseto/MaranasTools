from KBaseReport.KBaseReportClient import KBaseReport
from DataFileUtil.DataFileUtilClient import DataFileUtil
import uuid
import os
import errno
import shutil


class CreateReport(object):
    def __init__(self, callback_url, scratch_dir):
        self.scratch_dir = scratch_dir
        self.report_client = KBaseReport(callback_url)
        self.dfu_client = DataFileUtil(callback_url)

    def run(self, params):
        """
        params should have the following:
        * workspace_id = id of the workspace to save the report
        * text_input = some random string for the report
        * checkbox_input = a checkbox so we can add other params
        """
        ws_name = params.get('workspace_name', None)
        if ws_name is None:
            raise ValueError('workspace_name is required!')
        report_dir = self._build_report(params.get('text_input', ''),
                                        params.get('checkbox_input', 0))
        results = self._upload_report(report_dir, ws_name)
        return results

    def _build_report(self, text, checked):
        # make the directory
        res_dir = os.path.join(self.scratch_dir, "optstoic_out")
        report_dir = os.path.join(self.scratch_dir, str(uuid.uuid4()))
        self._mkdir_p(report_dir)
        # shutil.copy('/kb/module/lib/MaranasTools/Kbase_Logo_web.png', report_dir)
        shutil.copy(res_dir + "/pathway_001.png", report_dir)
        # make some example content
        report_file_content = """
        <html>
        <body>
        <div>
            I am an HTML report. Here's my contents:
        </div>
        <div>
            <b>Text input:</b>
            <pre>
            {}
            </pre>
        </div>
        <div>
            <b>Checkbox checked?</b>
            <br>
            {}
        </div>
        <div>
            <b>Here's a picture</b>
            <br>
            <img src="pathway_001.png" width="500px"/>
        </div>
        </body>
        </html>
        """.format(text, "yes" if checked == 1 else "no")

        # write the file to <report_dir>/index.html
        report_filename = os.path.join(report_dir, "index.html")
        with open(report_filename, "w") as report_file:
            report_file.write(report_file_content)
        return report_dir

    def _upload_report(self, report_dir, workspace_name):
        upload_info = self.dfu_client.file_to_shock({
            'file_path': report_dir,
            'pack': 'zip'
        })
        shock_id = upload_info['shock_id']

        report_params = {
            'message': 'This is a report message',
            'direct_html_link_index': 0,
            'html_links': [{
                'shock_id': shock_id,
                'name': 'index.html',
                'description': 'Just an example report'
            }],
            'report_object_name': 'NarrativeTest.example_report' + str(uuid.uuid4()),
            'workspace_name': workspace_name
        }
        report = self.report_client.create_extended_report(report_params)
        return {
            'report_ref': report['ref'],
            'report_name': report['name']
        }

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise