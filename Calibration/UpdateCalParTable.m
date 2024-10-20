function [CAL_table,showTable] = UpdateCalParTable(CAL_table,calType)

switch calType
    case 'UniGE TA'
        CAL_table.Data        = [{350};{740};{530};{1.1}];
        CAL_table.RowName     = ["Min WL (nm)"; "Max WL (nm)"; "Ctr. WL (nm)"; "pix/nm"];
        CAL_table.ColumnName  = "Det 1";
        showTable             = 1;
    case 'UniGE nsTA'
        CAL_table.Data        = [{350};{740};{530};{1.1}];
        CAL_table.RowName     = ["Min WL (nm)"; "Max WL (nm)"; "Ctr. WL (nm)"; "pix/nm"];
        CAL_table.ColumnName  = "Det 1";
        showTable             = 1;
    case 'UniGE TRUVIS-II (Intensity)'
        CAL_table.Data        = [{330};{720};{525};{1.25}];
        CAL_table.RowName     = ["Min WL (nm)"; "Max WL (nm)"; "Ctr. WL (nm)"; "pix/nm"];
        CAL_table.ColumnName  = "Det 1";
        showTable             = 1;
    case 'UoS TRIR'
        CAL_table.Data        = [{0,0};{0,0};{0,0};{0,0}];
        CAL_table.RowName     = ["Min cm-1"; "Max cm-1"; "Ctr. cm-1"; "pix/nm"];
        CAL_table.ColumnName  = ["Det 1"; "Det 2"];
        showTable             = 1;
    case {'RAL LIFEtime (Absorbance)','RAL LIFEtime (Intensity)'}
        CAL_table.Data        = [{0,0};{0,0};{0,0};{0,0}];
        CAL_table.RowName     = ["Min cm-1"; "Max cm-1"; "Ctr. cm-1"; "pix/nm"];
        CAL_table.ColumnName  = ["Det 1"; "Det 2"];
        showTable             = 1;
    otherwise
        showTable             = 0;
end

