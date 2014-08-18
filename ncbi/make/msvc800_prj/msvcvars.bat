@echo off
REM
REM $Id: msvcvars.bat,v 1.1 2006/09/29 20:30:06 ivanov Exp $
REM

@if not "%VSINSTALLDIR%"=="" goto devenv
@call "%VS80COMNTOOLS%vsvars32.bat"

:devenv

if exist "%VS80COMNTOOLS%..\IDE\VCExpress.*" set DEVENV="%VS80COMNTOOLS%..\IDE\VCExpress"
if exist "%VS80COMNTOOLS%..\IDE\devenv.*" set DEVENV="%VS80COMNTOOLS%..\IDE\devenv"

:end
