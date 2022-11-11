import csv
import os
import traceback
from threading import Thread
from pathlib import Path

import wx
from wx import xrc
import wx.grid
import wx.lib.newevent
import numpy as np

import config
import convert
import const
import utils
import view
import xtb
import xyzutils

CalcEndEvent, EVT_CALC_END = wx.lib.newevent.NewEvent()


class CalculationThread(Thread):

    def __init__(self, input_xyz_file: Path, xtb_params: xtb.XTBParams,
                 scans: list, constrains: list, force_constant, concerted: bool, keep_log: int,
                 num_threads: int, memory_per_thread: str, parent_window):
        Thread.__init__(self)
        self.input_xyz_file = input_xyz_file
        self.stop_file: Path = self.input_xyz_file.parent / (self.input_xyz_file.stem + const.STOP_FILE_SUFFIX)
        self.xtb_params = xtb_params
        self.scans = scans
        if len(scans) == 0:
            self.calculation_type = 'opt'
        else:
            self.calculation_type = 'scan'
        self.constrains = constrains
        self.force_constant = str(force_constant)
        self.concerted = concerted
        self.keep_log = keep_log
        self.num_threads = num_threads
        self.memory_per_thread = memory_per_thread
        self.parent_window = parent_window
        self.start()

    def run(self):
        if self.stop_file.exists():
            self.stop_file.unlink()

        try:
            xtb.setenv(self.num_threads, self.memory_per_thread)
            xtb.xtbscan(input_xyz_file=self.input_xyz_file, xtb_params=self.xtb_params,
                        scans=self.scans, constrains=self.constrains, force_constant=self.force_constant,
                        concerted=self.concerted, keep_log=self.keep_log)
        except:
            wx.PostEvent(self.parent_window, CalcEndEvent(input_xyz_file=self.input_xyz_file,
                                                          calculation_type=self.calculation_type,
                                                          success=False))
            raise
        else:
            wx.PostEvent(self.parent_window, CalcEndEvent(input_xyz_file=self.input_xyz_file,
                                                          calculation_type=self.calculation_type,
                                                          success=True))

    def terminate(self):
        self.stop_file.touch(exist_ok=True)


class MyFileDropTarget(wx.FileDropTarget):
    def __init__(self, window):
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, file_names):
        file = Path(file_names[0]).resolve()
        if file.suffix.lower() == '.csv':
            self.window.load_result_csv_file(file)
        else:
            self.window.load_input_xyz_file(file)
        return True


class CSVViewFrame(wx.Frame):
    def __init__(self, parent, title, data: list):
        wx.Frame.__init__(self, parent, -1, title)
        self.data = data
        self.init_frame()
        self.show_data()

    def init_frame(self):
        panel = wx.Panel(self, wx.ID_ANY)
        layout = wx.BoxSizer(wx.VERTICAL)
        self.grid_data = wx.grid.Grid(panel, wx.ID_ANY)
        font = wx.Font(12, wx.FONTFAMILY_MODERN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        self.grid_data.SetFont(font)
        layout.Add(self.grid_data, 1, wx.EXPAND | wx.ALL, border=3)
        panel.SetSizerAndFit(layout)

        self.grid_data.Bind(wx.EVT_KEY_DOWN, self.on_key_down_grid)

    def show_data(self):
        num_rows = len(self.data)
        num_cols = len(self.data[0])
        self.grid_data.CreateGrid(num_rows-1, num_cols-1)
        # Label for rows (number)
        for i in range(num_rows-1):
            self.grid_data.SetRowLabelValue(i, self.data[i+1][0])
        for i in range(num_cols-1):
            self.grid_data.SetColLabelValue(i, self.data[0][i+1])
        # show data
        for r in range(1, num_rows):
            for c in range(1, num_cols):
                self.grid_data.SetCellValue(r-1, c-1, self.data[r][c])
        self.grid_data.AutoSize()

        width = min(self.grid_data.Size[0] + 100, 1200)
        height = min(self.grid_data.Size[1] + 100, 800)
        self.SetSize((width, height))

    def copy_data(self):

        # Not selected but focused on somecell
        if len(self.grid_data.GetSelectionBlockBottomRight()) == 0:
            data = str(self.grid_data.GetCellValue(self.grid_data.GetGridCursorRow(),
                                                   self.grid_data.GetGridCursorCol()))

        else:
            rows = self.grid_data.GetSelectionBlockBottomRight()[0][0] - \
                   self.grid_data.GetSelectionBlockTopLeft()[0][0] + 1
            cols = self.grid_data.GetSelectionBlockBottomRight()[0][1] - \
                   self.grid_data.GetSelectionBlockTopLeft()[0][1] + 1

            data = ''
            for r in range(rows):
                for c in range(cols):
                    data = data + str(self.grid_data.GetCellValue(self.grid_data.GetSelectionBlockTopLeft()[0][0] + r,
                                                                  self.grid_data.GetSelectionBlockTopLeft()[0][1] + c))
                    if c < cols - 1:
                        data = data + '\t'
                data = data + '\n'

        clipboard = wx.TextDataObject()
        clipboard.SetText(data)
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            pass

    def on_key_down_grid(self, event):
        if event.ControlDown() and event.GetKeyCode() == 67:
            self.copy_data()


class XTBScanApp(wx.App):

    # initialization
    def OnInit(self):
        self.resource = xrc.XmlResource('./xtbscan.xrc')
        self.init_frame()

        # internal variables
        self.current_input_xyz_file: Path = None
        self.current_input_xyz_coordinates: np.ndarray = None
        self.current_result_csv_file: Path = None
        self.current_result_xyz_file: Path = None
        self.current_result_xyz_atoms: np.ndarray = None
        self.current_result_xyz_coordinates: list = None

        # calculation thread: None means vacant
        self.calculation_thread = None

        return True

    def init_frame(self):
        self.frame = self.resource.LoadFrame(None, 'frame')
        self.get_controls_from_xrc()
        self.set_events()
        self.init_controls()

        # Drug & Drop settings
        dt = MyFileDropTarget(self)
        self.frame.SetDropTarget(dt)

        self.frame.Show()

    def get_controls_from_xrc(self):
        # Codes for loading objects
        self.text_ctrl_input_xyz_file: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_input_xyz_file')
        self.button_input_xyz_file: wx.Button = xrc.XRCCTRL(self.frame, 'button_input_xyz_file')
        self.button_view_input_xyz_file: wx.Button = xrc.XRCCTRL(self.frame, 'button_view_input_xyz_file')
        self.text_ctrl_input_xyz_coordinates: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_input_xyz_coordinates')
        self.text_ctrl_result_csv_file: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_result_csv_file')
        self.button_result_csv_file: wx.Button = xrc.XRCCTRL(self.frame, 'button_result_csv_file')
        self.button_view_result_xyz: wx.Button = xrc.XRCCTRL(self.frame, 'button_view_result_xyz')
        self.button_view_result_csv: wx.Button = xrc.XRCCTRL(self.frame, 'button_view_result_csv')
        self.text_ctrl_view_result_number: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_view_result_number')
        self.button_view_xyz_selected: wx.Button = xrc.XRCCTRL(self.frame, 'button_view_xyz_selected')
        self.button_view_copy_selected: wx.Button = xrc.XRCCTRL(self.frame, 'button_view_copy_selected')
        self.button_plot_result: wx.Button = xrc.XRCCTRL(self.frame, 'button_plot_result')
        self.checkbox_plot_annotation: wx.CheckBox = xrc.XRCCTRL(self.frame, 'checkbox_plot_annotation')
        self.button_plot_surface: wx.Button = xrc.XRCCTRL(self.frame, 'button_plot_surface')
        self.choice_xtb_method: wx.Choice = xrc.XRCCTRL(self.frame, 'choice_xtb_method')
        self.text_ctrl_xtb_charge: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_xtb_charge')
        self.text_ctrl_xtb_uhf: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_xtb_uhf')
        self.radio_box_xtb_solvation: wx.RadioBox = xrc.XRCCTRL(self.frame, 'radio_box_xtb_solvation')
        self.choice_xtb_solvent: wx.Choice = xrc.XRCCTRL(self.frame, 'choice_xtb_solvent')
        self.text_ctrl_force_constant: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_force_constant')
        self.text_ctrl_scan_cpus: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_cpus')
        self.text_ctrl_scan_memory: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_memory')
        self.choice_keep: wx.Choice = xrc.XRCCTRL(self.frame, 'choice_keep')
        self.button_run: wx.Button = xrc.XRCCTRL(self.frame, 'button_run')
        self.checkbox_scan_concerted: wx.CheckBox = xrc.XRCCTRL(self.frame, 'checkbox_scan_concerted')
        self.list_ctrl_scans: wx.ListCtrl = xrc.XRCCTRL(self.frame, 'list_ctrl_scans')
        self.button_scan_remove: wx.Button = xrc.XRCCTRL(self.frame, 'button_scan_remove')
        self.button_scan_add: wx.Button = xrc.XRCCTRL(self.frame, 'button_scan_add')
        self.choice_scan_type: wx.Choice = xrc.XRCCTRL(self.frame, 'choice_scan_type')
        self.text_ctrl_scan_atoms: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_atoms')
        self.text_ctrl_scan_start: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_start')
        self.text_ctrl_scan_end: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_end')
        self.text_ctrl_scan_num_step: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_scan_num_step')
        self.button_scan_get_current: wx.Button = xrc.XRCCTRL(self.frame, 'button_scan_get_current')
        self.list_ctrl_constrains: wx.ListCtrl = xrc.XRCCTRL(self.frame, 'list_ctrl_constrains')
        self.button_constrain_remove: wx.Button = xrc.XRCCTRL(self.frame, 'button_constrain_remove')
        self.button_constrain_add: wx.Button = xrc.XRCCTRL(self.frame, 'button_constrain_add')
        self.choice_constrain_type: wx.Choice = xrc.XRCCTRL(self.frame, 'choice_constrain_type')
        self.text_ctrl_constrain_atoms: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_constrain_atoms')
        self.text_ctrl_constrain_value: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_constrain_value')
        self.button_constrain_get_current: wx.Button = xrc.XRCCTRL(self.frame, 'button_constrain_get_current')
        self.text_ctrl_log: wx.TextCtrl = xrc.XRCCTRL(self.frame, 'text_ctrl_log')

        # Codes for checking objects
        assert self.text_ctrl_input_xyz_file is not None
        assert self.button_input_xyz_file is not None
        assert self.button_view_input_xyz_file is not None
        assert self.text_ctrl_input_xyz_coordinates is not None
        assert self.text_ctrl_result_csv_file is not None
        assert self.button_result_csv_file is not None
        assert self.button_view_result_xyz is not None
        assert self.button_view_result_csv is not None
        assert self.text_ctrl_view_result_number is not None
        assert self.button_view_xyz_selected is not None
        assert self.button_view_copy_selected is not None
        assert self.button_plot_result is not None
        assert self.checkbox_plot_annotation is not None
        assert self.button_plot_surface is not None
        assert self.choice_xtb_method is not None
        assert self.text_ctrl_xtb_charge is not None
        assert self.text_ctrl_xtb_uhf is not None
        assert self.radio_box_xtb_solvation is not None
        assert self.choice_xtb_solvent is not None
        assert self.text_ctrl_force_constant is not None
        assert self.text_ctrl_scan_cpus is not None
        assert self.text_ctrl_scan_memory is not None
        assert self.choice_keep is not None
        assert self.button_run is not None
        assert self.checkbox_scan_concerted is not None
        assert self.list_ctrl_scans is not None
        assert self.button_scan_remove is not None
        assert self.button_scan_add is not None
        assert self.choice_scan_type is not None
        assert self.text_ctrl_scan_atoms is not None
        assert self.text_ctrl_scan_start is not None
        assert self.text_ctrl_scan_end is not None
        assert self.text_ctrl_scan_num_step is not None
        assert self.button_scan_get_current is not None
        assert self.list_ctrl_constrains is not None
        assert self.button_constrain_remove is not None
        assert self.button_constrain_add is not None
        assert self.choice_constrain_type is not None
        assert self.text_ctrl_constrain_atoms is not None
        assert self.text_ctrl_constrain_value is not None
        assert self.button_constrain_get_current is not None
        assert self.text_ctrl_log is not None

    def set_events(self):

        # Input file
        self.button_input_xyz_file.Bind(wx.EVT_BUTTON, self.on_button_input_xyz_file)
        self.text_ctrl_input_xyz_file.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter_input_xyz_file)

        # Result file
        self.button_result_csv_file.Bind(wx.EVT_BUTTON, self.on_button_result_csv_file)
        self.button_view_input_xyz_file.Bind(wx.EVT_BUTTON, self.on_button_view_input_xyz_file)
        self.text_ctrl_result_csv_file.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter_result_csv_file)
        self.button_view_result_xyz.Bind(wx.EVT_BUTTON, self.on_button_view_result_xyz)
        self.button_view_result_csv.Bind(wx.EVT_BUTTON, self.on_button_view_result_csv)
        self.button_view_xyz_selected.Bind(wx.EVT_BUTTON, self.on_button_view_xyz_selected)
        self.button_view_copy_selected.Bind(wx.EVT_BUTTON, self.on_button_view_copy_selected)
        self.button_plot_result.Bind(wx.EVT_BUTTON, self.on_button_plot_result)
        self.button_plot_surface.Bind(wx.EVT_BUTTON, self.on_button_plot_surface)

        # Scans
        self.button_scan_get_current.Bind(wx.EVT_BUTTON, self.on_button_scan_get_current)
        self.button_scan_add.Bind(wx.EVT_BUTTON, self.on_button_scan_add)
        self.button_scan_remove.Bind(wx.EVT_BUTTON, self.on_button_scan_remove)

        # Constrains
        self.button_constrain_get_current.Bind(wx.EVT_BUTTON, self.on_button_constrain_get_current)
        self.button_constrain_add.Bind(wx.EVT_BUTTON, self.on_button_constain_add)
        self.button_constrain_remove.Bind(wx.EVT_BUTTON, self.on_button_constain_remove)

        # Run
        self.button_run.Bind(wx.EVT_BUTTON, self.on_button_run)

        # Event when calculation finished
        self.Bind(EVT_CALC_END, self.on_calc_end)

    def init_controls(self):
        # XTB solvent
        for solvent in const.XTB_SOLVENT_LIST:
            self.choice_xtb_solvent.Append(solvent)

        # keep logs
        for keep in ['No', 'Yes', 'when fail']:
            self.choice_keep.Append(keep)
            self.choice_keep.SetStringSelection('when fail')

        # Scans
        for header in ['Type', 'Atoms', 'Start', 'End', 'Num Step']:
            self.list_ctrl_scans.AppendColumn(header)

        # Constrains
        for header in ['Type', 'Atoms', 'Value']:
            self.list_ctrl_constrains.AppendColumn(header)

    def with_logging_exceptions(func):
        def inner(self, *args, **kwargs):
            try:
                result = func(self, *args, **kwargs)
            except Exception as e:
                if config.DEBUG:
                    tb = list(traceback.TracebackException.from_exception(e).format())
                    self.logging(tb)
                else:
                    self.logging(e.args)
                return None
            else:
                return result

        return inner

    # main working functions (called by event-handlers)
    @with_logging_exceptions
    def load_input_xyz_file(self, file: Path):
        # check file type by suffix and convert Gaussian file to xyz
        if file.suffix.lower() in ['.gjf', '.gjc', '.com']:
            file = convert.gaussian_input_to_xyz(file)
            self.logging('Convert Gaussian input to xyz: {}'.format(file))
        elif file.suffix.lower() in ['.log', '.out']:
            file = convert.gaussian_log_to_xyz(file)
            self.logging('Convert Gaussian log (final structure) to xyz: {}'.format(file))

        self.text_ctrl_input_xyz_file.SetValue(str(file.resolve()))
        atoms, coordinates = xyzutils.read_single_xyz_file(file)
        self.text_ctrl_input_xyz_coordinates.SetValue(xyzutils.get_xyz_string(atoms, coordinates))
        self.current_input_xyz_coordinates = coordinates
        self.current_input_xyz_file = file.resolve()

    @with_logging_exceptions
    def load_result_csv_file(self, file: Path):
        self.text_ctrl_result_csv_file.SetValue(str(file.resolve()))
        self.current_result_csv_file = file.resolve()

        # check and read the corresponding xyz file
        xyz_file = file.resolve().parent / (file.stem + '.xyz')
        if not xyz_file.exists():
            self.logging('The related result xyz file not found: {:}'.format(xyz_file))
            self.current_result_xyz_file = None
            self.current_result_xyz_coordinates = None
            self.current_result_xyz_atoms = None
        else:
            self.logging('The related result xyz file is read: {:}'.format(xyz_file))
            self.load_result_xyz_file(xyz_file)

    def load_result_xyz_file(self, file: Path):
        self.current_result_xyz_file = file.resolve()
        self.current_result_xyz_atoms, self.current_result_xyz_coordinates, _ = xyzutils.read_sequential_xyz_file(file)

    def logging(self, message):
        """
        output logs
        """
        log_string = (''.join(message)).rstrip()
        self.text_ctrl_log.write(log_string + '\n')
        if config.DEBUG:
            print(log_string)

    # Event-handlers
    # Input event-handlers
    @with_logging_exceptions
    def on_button_input_xyz_file(self, event):
        dialog = wx.FileDialog(None, 'Select input structure file',
                               wildcard='XYZ or Gaussian file  (*.xyz;*.gjf;*.gjc;*.com;*.log;*.out)|'
                                        '*.xyz;*.gjf;*.gjc;*.com;*.log;*.out'
                                        '|All files (*.*)|*.*',
                               style=wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            file = Path(dialog.GetPath())
            self.load_input_xyz_file(file)
        dialog.Destroy()

    @with_logging_exceptions
    def on_text_enter_input_xyz_file(self, event):
        file = Path(self.text_ctrl_input_xyz_file.GetValue()).resolve()
        if file.exists():
            self.load_input_xyz_file(file)
        else:
            self.logging('Input file not found: {}'.format(file))
            self.current_input_xyz_file = None
            self.current_input_xyz_coordinates = None
            self.text_ctrl_input_xyz_coordinates.SetValue('')

    @with_logging_exceptions
    def on_button_view_input_xyz_file(self, event):
        if self.current_input_xyz_file is not None:
            view.view_xyz_file(self.current_input_xyz_file)
        else:
            self.logging('No input file')

    # Result event-handlers
    @with_logging_exceptions
    def on_button_result_csv_file(self, event):
        dialog = wx.FileDialog(None, 'Select result csv file',
                               wildcard='CSV file (*.csv)|*.csv|All files (*.*)|*.*',
                               style=wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            file = Path(dialog.GetPath())
            self.load_result_csv_file(file)
        dialog.Destroy()

    @with_logging_exceptions
    def on_text_enter_result_csv_file(self, event):
        file = Path(self.text_ctrl_result_csv_file.GetValue()).resolve()
        if file.exists():
            self.load_result_csv_file(file)
        else:
            self.logging('Input file not found: {}'.format(file))
            self.current_result_csv_file = None
            self.current_result_xyz_file = None
            self.current_result_xyz_coordinates = None
            self.current_result_xyz_atoms = None

    @with_logging_exceptions
    def on_button_view_result_xyz(self, event):
        if self.current_result_xyz_file is not None:
            view.view_xyz_file(self.current_result_xyz_file)
        else:
            self.logging('No result xyz file.')

    @with_logging_exceptions
    def on_button_view_result_csv(self, event):
        if self.current_result_csv_file is not None:

            with self.current_result_csv_file.open(mode='r') as f:
                reader = csv.reader(f)
                csv_data = list(reader)[1:]

            title = str(self.current_result_csv_file)

            # check: all lines have same length
            length = len(csv_data[0])
            for line in csv_data[1:]:
                if len(line) != length:
                    self.logging('CSV data is no valid. All data lines must have the same length.')
                    return
            # show
            csv_view_frame = CSVViewFrame(self.frame, title=title, data=csv_data)
            csv_view_frame.Show(True)

        else:
            self.logging('No result csv file.')

    @with_logging_exceptions
    def on_button_view_xyz_selected(self, event):
        if self.current_result_xyz_coordinates is None or self.current_result_xyz_atoms is None:
            self.logging('No valid xyz file')
            return

        # check input number
        try:
            num = int(self.text_ctrl_view_result_number.GetValue())
        except:
            self.logging('The given number is not valid.')
            return

        # check input number range
        num_structure = len(self.current_result_xyz_coordinates)
        if num <= 0 or num > num_structure:
            self.logging('The given number is out of range. Input 0 to {:}.'.format(num_structure))
            return

        view.view_xyz_structure(atoms=self.current_result_xyz_atoms,
                                coordinates=self.current_result_xyz_coordinates[num - 1],
                                title='#. {:}'.format(num))

    @with_logging_exceptions
    def on_button_view_copy_selected(self, event):
        if self.current_result_xyz_coordinates is None or self.current_result_xyz_atoms is None:
            self.logging('No valid xyz file.')
            return

        # check input number
        try:
            num = int(self.text_ctrl_view_result_number.GetValue())
        except:
            self.logging('The given number is not valid.')
            return

        # check input number range
        num_structure = len(self.current_result_xyz_coordinates)
        if num <= 0 or num > num_structure:
            self.logging('The given number is out of range. Input 1 to {:}.'.format(num_structure))
            return

        xyz_string = xyzutils.get_xyz_string(self.current_result_xyz_atoms,
                                             self.current_result_xyz_coordinates[num - 1])

        clipboard = wx.TextDataObject()
        clipboard.SetText(xyz_string)
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
            self.logging('XYZ coordinates copied to clipboard.')
        else:
            self.logging('Clipboard is not accessible.')


    @with_logging_exceptions
    def on_button_plot_result(self, event):
        if self.current_result_csv_file is not None:
            annotation = self.checkbox_plot_annotation.GetValue()
            view.plot_scan(self.current_result_csv_file, annotation=annotation)
        else:
            self.logging('No result csv file.')

    @with_logging_exceptions
    def on_button_plot_surface(self, event):
        if self.current_result_csv_file is not None:
            view.plot_surface(self.current_result_csv_file)
        else:
            self.logging('No result csv file.')

    @with_logging_exceptions
    def on_button_scan_get_current(self, event):
        if self.current_input_xyz_coordinates is None:
            self.logging('No input coordinates.')
            return
        stype = self.choice_scan_type.GetStringSelection()
        atoms = utils.read_atom_list_string(self.text_ctrl_scan_atoms.GetValue())

        # check atoms in coordinates
        for i in atoms:
            if i >= len(self.current_input_xyz_coordinates):
                self.logging('Input cooridinates only have {} atoms.'.format(len(self.current_input_xyz_coordinates)))
                return

        if stype == 'distance':
            if len(atoms) != 2:
                self.logging('2 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_distance(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_scan_start.SetValue(v)

        if stype == 'angle':
            if len(atoms) != 3:
                self.logging('3 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_angle(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_scan_start.SetValue(v)

        if stype == 'dihedral':
            if len(atoms) != 4:
                self.logging('4 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_dihedral(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_scan_start.SetValue(v)

    @with_logging_exceptions
    def on_button_constrain_get_current(self, event):
        if self.current_input_xyz_coordinates is None:
            self.logging('No input coordinates.')
            return
        ctype = self.choice_constrain_type.GetStringSelection()
        atoms = utils.read_atom_list_string(self.text_ctrl_constrain_atoms.GetValue())

        # check atoms in coordinates
        for i in atoms:
            if i >= len(self.current_input_xyz_coordinates):
                self.logging('Input cooridinates only have {} atoms.'.format(len(self.current_input_xyz_coordinates)))
                return

        if ctype == 'distance':
            if len(atoms) != 2:
                self.logging('2 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_distance(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_constrain_value.SetValue(v)

        if ctype == 'angle':
            if len(atoms) != 3:
                self.logging('3 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_angle(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_constrain_value.SetValue(v)

        if ctype == 'dihedral':
            if len(atoms) != 4:
                self.logging('4 atoms should be given to calculate distance.')
                return
            else:
                v = xyzutils.calc_dihedral(self.current_input_xyz_coordinates, atoms, string=True)
                self.text_ctrl_constrain_value.SetValue(v)

    @with_logging_exceptions
    def on_button_scan_add(self, event):
        stype = self.choice_scan_type.GetStringSelection()
        atoms = self.text_ctrl_scan_atoms.GetValue().strip()
        start = self.text_ctrl_scan_start.GetValue().strip()
        end = self.text_ctrl_scan_end.GetValue().strip()
        num_step = self.text_ctrl_scan_num_step.GetValue().strip()
        errors = scan_validation(stype, atoms, start, end, num_step, None)

        if len(errors) == 0:
            atoms = utils.expand_atom_list_string(atoms)
            self.list_ctrl_scans.Append([stype, atoms, start, end, num_step])
            self.text_ctrl_scan_atoms.SetValue('')
            self.text_ctrl_scan_start.SetValue('')
            self.text_ctrl_scan_end.SetValue('')
            self.text_ctrl_scan_num_step.SetValue('')
        else:
            self.logging('\n'.join(errors))

    @with_logging_exceptions
    def on_button_scan_remove(self, event):
        for i in reversed(range(self.list_ctrl_scans.GetItemCount())):
            if self.list_ctrl_scans.IsSelected(i):
                self.list_ctrl_scans.DeleteItem(i)

    @with_logging_exceptions
    def on_button_constain_add(self, event):
        ctype = self.choice_constrain_type.GetStringSelection()
        atoms = self.text_ctrl_constrain_atoms.GetValue().strip()
        value = self.text_ctrl_constrain_value.GetValue().lower().strip()
        errors = constrain_validation(ctype, atoms, value, None)

        if len(errors) == 0:
            atoms = utils.expand_atom_list_string(atoms)
            if value == '':
                value = 'auto'
            self.list_ctrl_constrains.Append([ctype, atoms, value])
            self.text_ctrl_constrain_atoms.SetValue('')
            self.text_ctrl_constrain_value.SetValue('auto')
        else:
            self.logging('\n'.join(errors))

    @with_logging_exceptions
    def on_button_constain_remove(self, event):
        for i in reversed(range(self.list_ctrl_constrains.GetItemCount())):
            if self.list_ctrl_constrains.IsSelected(i):
                self.list_ctrl_constrains.DeleteItem(i)

    @with_logging_exceptions
    def on_button_run(self, event):
        if self.calculation_thread is not None:

            # Termination check
            msgbox = wx.MessageDialog(None,
                                      'Stop the running calculation?',
                                      'Calculation in running',
                                      style=wx.YES_NO)
            stop_check = msgbox.ShowModal()
            if stop_check == wx.ID_YES:
                msgbox.Destroy()
                self.calculation_thread.terminate()
                self.logging('Termination requested.')
            else:
                msgbox.Destroy()
                return

        else:
            # Input xyz
            if self.current_input_xyz_file is None or self.current_input_xyz_coordinates is None:
                self.logging('No input file.')
                return

            # Settings for calulcation
            # XTB Parameters
            xtb_method = self.choice_xtb_method.GetStringSelection()
            try:
                xtb_charge = int(self.text_ctrl_xtb_charge.GetValue().strip())
            except ValueError:
                self.logging('Charge is not valid.')
                return
            try:
                xtb_uhf = int(self.text_ctrl_xtb_uhf.GetValue().strip())
                assert xtb_uhf >= 0
            except:
                self.logging('UHF is not valid.')
                return
            xtb_solvation = self.radio_box_xtb_solvation.GetStringSelection().lower()
            if xtb_solvation == 'none':
                xtb_solvation = None
                xtb_solvent = ''
            else:
                xtb_solvent = self.choice_xtb_solvent.GetStringSelection().strip().lower()
                if not xtb_solvent:
                    self.logging('Solvent must be selected.')
                    return
            xtb_params = xtb.XTBParams(method=xtb_method, charge=xtb_charge, uhf=xtb_uhf,
                                       solvation=xtb_solvation, solvent=xtb_solvent)

            force_constant = self.text_ctrl_force_constant.GetValue().strip()
            try:
                _fc = float(force_constant)
                assert _fc > 0.0
            except:
                self.logging('Force Constant is not valid. It should positive value (0.5, 1.0 etc.)')
                return

            # cpus, memory, keep_log
            try:
                num_cpus = int(self.text_ctrl_scan_cpus.GetValue().strip())
                assert num_cpus >= 1
            except:
                self.logging('The Number of CPUs is not valid. It must be positive integer.')
                return
            memory_per_cpu: str = self.text_ctrl_scan_memory.GetValue().strip()
            memory_per_cpu = memory_per_cpu.upper().replace(' ', '').replace('　', '')
            keep_log = self.choice_keep.GetStringSelection().lower()
            if keep_log == 'yes':
                keep_log = 2
            elif keep_log == 'when fail':
                keep_log = 1
            else:
                keep_log = 0

            # Scans
            scans = []
            for n in range(self.list_ctrl_scans.GetItemCount()):
                stype = self.list_ctrl_scans.GetItem(n, 0).GetText()
                atoms = self.list_ctrl_scans.GetItem(n, 1).GetText()
                start = self.list_ctrl_scans.GetItem(n, 2).GetText()
                end = self.list_ctrl_scans.GetItem(n, 3).GetText()
                num_step = self.list_ctrl_scans.GetItem(n, 4).GetText()
                errors = scan_validation(stype, atoms, start, end, num_step, self.current_input_xyz_coordinates)
                if len(errors) > 0:
                    self.logging('Errors in Scan line {:}'.format(n + 1))
                    self.logging('\n'.join(errors))
                    return
                else:
                    atom_indices = utils.read_atom_list_string(atoms)
                    scans.append(xtb.XTBScan(stype, atom_indices, start, end, int(num_step)))

            # Constrains
            constrains = []
            for n in range(self.list_ctrl_constrains.GetItemCount()):
                ctype = self.list_ctrl_constrains.GetItem(n, 0).GetText()
                atoms = self.list_ctrl_constrains.GetItem(n, 1).GetText()
                value = self.list_ctrl_constrains.GetItem(n, 2).GetText()
                errors = constrain_validation(ctype, atoms, value, self.current_input_xyz_coordinates)
                if len(errors) > 0:
                    self.logging('Errors in Constrain line {:}'.format(n + 1))
                    self.logging('\n'.join(errors))
                    return
                else:
                    atom_indices = utils.read_atom_list_string(atoms)
                    constrains.append(xtb.XTBConstrain(ctype, atom_indices, value))

            # Start Calculation
            self.calculation_thread = CalculationThread(input_xyz_file=self.current_input_xyz_file,
                                                        xtb_params=xtb_params,
                                                        scans=scans,
                                                        constrains=constrains,
                                                        force_constant=force_constant,
                                                        concerted=self.checkbox_scan_concerted.IsChecked(),
                                                        keep_log=keep_log,
                                                        num_threads=num_cpus,
                                                        memory_per_thread=memory_per_cpu,
                                                        parent_window=self)

            self.button_run.SetLabelText('Stop')
            self.logging('Calculation started: {}'.format(self.current_input_xyz_file))

    @with_logging_exceptions
    def on_calc_end(self, event: CalcEndEvent):

        self.calculation_thread = None
        self.button_run.SetLabelText('Run')

        input_xyz_file: Path = event.input_xyz_file

        if event.success:
            self.logging('Calculation successfully finished: {}'.format(input_xyz_file))
            if event.calculation_type == 'opt':
                pass
            else:
                msgbox = wx.MessageDialog(None,
                                          'Calculation successfully finished: {}\n'.format(input_xyz_file.name) + \
                                          'Load result files?',
                                          'Calculation Finished',
                                          style=wx.YES_NO)
                load_check = msgbox.ShowModal()
                if load_check == wx.ID_YES:
                    msgbox.Destroy()
                    result_csv_file = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_CSV_SUFFIX)
                    self.load_result_csv_file(result_csv_file)
                else:
                    msgbox.Destroy()
        else:
            self.logging('Calculation failed or terminated: {}'.format(input_xyz_file))


def scan_validation(stype: str, atoms: str, start: str, end: str, num_step: str, coodrinates=None):
    """
    if input are not valid, return error messages in list [str].
    If ok, return empty list [].
    """
    errors = []

    if stype.lower() not in ['distance', 'angle', 'dihedral']:
        errors.append('Type must be distance, angle, or dihedral.')

    try:
        atom_indices = utils.read_atom_list_string(atoms)
        if stype == 'distance':
            if len(atom_indices) != 2:
                errors.append('2 atoms should be given for distance.')
        if stype == 'angle':
            if len(atom_indices) != 3:
                errors.append('3 atoms should be given for angle.')
        if stype == 'dihedral':
            if len(atom_indices) != 4:
                errors.append('4 atoms should be given for dihedral.')
        if coodrinates is not None:
            for i in atom_indices:
                if i >= len(coodrinates):
                    errors.append('Given atom index is out of range.')
                    break
    except:
        errors.append('Atoms are not valid.')

    try:
        float(start)
    except ValueError:
        errors.append('Start is not valid value.')

    try:
        float(end)
    except ValueError:
        errors.append('End is not valid value.')

    try:
        num_step = int(num_step)
    except ValueError:
        errors.append('Num step is not valid value. It must be integer')
    else:
        if num_step < 2:
            errors.append('Num step must be >= 2.')
            return

    return errors


def constrain_validation(ctype: str, atoms: str, value: str, coodrinates=None):
    """
    if input are not valid, return error messages in list [str].
    If ok, return empty list [].
    """
    errors = []

    if ctype.lower() not in ['distance', 'angle', 'dihedral']:
        errors.append('Type must be distance, angle, or dihedral.')

    try:
        atom_indices = utils.read_atom_list_string(atoms)
        if ctype == 'distance':
            if len(atom_indices) != 2:
                errors.append('2 atoms should be given for distance.')
        if ctype == 'angle':
            if len(atom_indices) != 3:
                errors.append('3 atoms should be given for angle.')
        if ctype == 'dihedral':
            if len(atom_indices) != 4:
                errors.append('4 atoms should be given for dihedral.')
        if coodrinates is not None:
            for i in atom_indices:
                if i >= len(coodrinates):
                    errors.append('Given atom index is out of range.')
                    break
    except:
        errors.append('Atoms are not valid.')

    try:
        if value.lower() not in ['', 'auto']:
            float(value)
    except ValueError:
        errors.append('Constrain value is not valid.')

    return errors


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    app = XTBScanApp(False)
    app.MainLoop()
