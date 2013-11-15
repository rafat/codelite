##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=qrdecomp
ConfigurationName      :=Debug
WorkspacePath          := "C:\Documents and Settings\USER\My Documents\Downloads\codelite\ctest"
ProjectPath            := "C:\Documents and Settings\USER\My Documents\Downloads\codelite\ctest\qrdecomp"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=USER
Date                   :=9/3/2013
CodeLitePath           :="C:\Program Files\CodeLite"
LinkerName             :=gcc
SharedObjectLinkerName :=gcc -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="qrdecomp.txt"
PCHCompileFlags        :=
MakeDirCommand         :=makedir
RcCmpOptions           := 
RcCompilerName         :=windres
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := gcc
CC       := gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)


##
## User defined environment variables
##
CodeLiteDir:=C:\Program Files\CodeLite
Objects0=$(IntermediateDirectory)/main$(ObjectSuffix) $(IntermediateDirectory)/qrdecomp$(ObjectSuffix) $(IntermediateDirectory)/matrix$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

$(IntermediateDirectory)/.d:
	@$(MakeDirCommand) "./Debug"

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main$(ObjectSuffix): main.c $(IntermediateDirectory)/main$(DependSuffix)
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/qrdecomp/main.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main$(DependSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main$(ObjectSuffix) -MF$(IntermediateDirectory)/main$(DependSuffix) -MM "main.c"

$(IntermediateDirectory)/main$(PreprocessSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main$(PreprocessSuffix) "main.c"

$(IntermediateDirectory)/qrdecomp$(ObjectSuffix): qrdecomp.c $(IntermediateDirectory)/qrdecomp$(DependSuffix)
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/qrdecomp/qrdecomp.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/qrdecomp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/qrdecomp$(DependSuffix): qrdecomp.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/qrdecomp$(ObjectSuffix) -MF$(IntermediateDirectory)/qrdecomp$(DependSuffix) -MM "qrdecomp.c"

$(IntermediateDirectory)/qrdecomp$(PreprocessSuffix): qrdecomp.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/qrdecomp$(PreprocessSuffix) "qrdecomp.c"

$(IntermediateDirectory)/matrix$(ObjectSuffix): matrix.c $(IntermediateDirectory)/matrix$(DependSuffix)
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/qrdecomp/matrix.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/matrix$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/matrix$(DependSuffix): matrix.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/matrix$(ObjectSuffix) -MF$(IntermediateDirectory)/matrix$(DependSuffix) -MM "matrix.c"

$(IntermediateDirectory)/matrix$(PreprocessSuffix): matrix.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/matrix$(PreprocessSuffix) "matrix.c"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) $(IntermediateDirectory)/main$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/main$(DependSuffix)
	$(RM) $(IntermediateDirectory)/main$(PreprocessSuffix)
	$(RM) $(IntermediateDirectory)/qrdecomp$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/qrdecomp$(DependSuffix)
	$(RM) $(IntermediateDirectory)/qrdecomp$(PreprocessSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(DependSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(PreprocessSuffix)
	$(RM) $(OutputFile)
	$(RM) $(OutputFile).exe
	$(RM) "../.build-debug/qrdecomp"


