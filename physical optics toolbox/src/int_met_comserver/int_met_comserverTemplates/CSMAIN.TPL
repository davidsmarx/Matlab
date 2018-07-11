@#nooutput
@#library "AGCallouts.dll"
@(TYPELIB_ID := GUIDGEN()@)
@#template "server.idl"
@#template "serverglobal.f90|@GLOBAL.NAME@Global.f90"
@#template "serverhelper.f90"
@#IFTRUE EXESERVER
@#template "exemain.f90"
@#ELSE
@#template "server.def"
@#template "dllmain.f90"
@#ENDIF
@#template "clsfactty.f90"
@#template "clsfact.f90"
@#per CLASS
@#template "ccccty.f90|@CLASS.SHORTNAME@TY.f90"
@#template "userty.f90|U@CLASS.SHORTNAME@TY.f90|U@CLASS.SHORTNAME@TY+.f90"
@#endper
