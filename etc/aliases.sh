#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     aliases.sh
#
# Description
#     Aliases for working with OpenFOAM
#     Sourced from OpenFOAM-<VERSION>/etc/bashrc and/or ~/.bashrc
#
#------------------------------------------------------------------------------

# Change compiled version aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wmSET='. $WM_PROJECT_DIR/etc/bashrc'
alias wm64='wmSET WM_ARCH_OPTION=64'
alias wm32='wmSET WM_ARCH_OPTION=32'
alias wmSP='wmSET WM_PRECISION_OPTION=SP'
alias wmDP='wmSET WM_PRECISION_OPTION=DP'

# Toggle wmakeScheduler on/off
#  - also need to set WM_HOSTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wmSchedON='export WM_SCHEDULER=$WM_PROJECT_DIR/wmake/wmakeScheduler'
alias wmSchedOFF='unset WM_SCHEDULER'

# Change paraview version
# ~~~~~~~~~~~~~~~~~~~~~~~
unset foamPV
foamPV()
{
    export ParaView_VERSION=$1
    . $WM_PROJECT_DIR/etc/apps/paraview3/bashrc
    echo "paraview-$ParaView_VERSION  (major: $ParaView_MAJOR)"
    echo "dir: $ParaView_DIR"
    [ -d "$ParaView_DIR" ] || echo "WARNING: directory does not exist"
}


# Change directory aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~
alias src='cd $FOAM_SRC'
alias lib='cd $FOAM_LIB'
alias run='cd $FOAM_RUN'
alias foam='cd $WM_PROJECT_DIR'
alias foamsrc='cd $FOAM_SRC/$WM_PROJECT'
alias foamfv='cd $FOAM_SRC/finiteVolume'
alias app='cd $FOAM_APP'
alias util='cd $FOAM_UTILITIES'
alias sol='cd $FOAM_SOLVERS'
alias tut='cd $FOAM_TUTORIALS'

alias foamApps='cd $FOAM_APP'
alias foamSol='cd $FOAM_SOLVERS'
alias foamTuts='cd $FOAM_TUTORIALS'
alias foamUtils='cd $FOAM_UTILITIES'
alias foam3rdParty='cd $WM_THIRD_PARTY_DIR'

# -----------------------------------------------------------------------------
