##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=newtonmin
ConfigurationName      :=Debug
WorkspacePath          := "C:\Documents and Settings\USER\My Documents\Downloads\codelite\ctest"
ProjectPath            := "C:\Documents and Settings\USER\My Documents\Downloads\codelite\ctest\newtonmin"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=USER
Date                   :=11/8/2013
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
ObjectsFileList        :="newtonmin.txt"
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
Objects0=$(IntermediateDirectory)/main$(ObjectSuffix) $(IntermediateDirectory)/matrix$(ObjectSuffix) $(IntermediateDirectory)/newtonmin$(ObjectSuffix) 



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
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/newtonmin/main.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main$(DependSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main$(ObjectSuffix) -MF$(IntermediateDirectory)/main$(DependSuffix) -MM "main.c"

$(IntermediateDirectory)/main$(PreprocessSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main$(PreprocessSuffix) "main.c"

$(IntermediateDirectory)/matrix$(ObjectSuffix): matrix.c $(IntermediateDirectory)/matrix$(DependSuffix)
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/newtonmin/matrix.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/matrix$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/matrix$(DependSuffix): matrix.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/matrix$(ObjectSuffix) -MF$(IntermediateDirectory)/matrix$(DependSuffix) -MM "matrix.c"

$(IntermediateDirectory)/matrix$(PreprocessSuffix): matrix.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/matrix$(PreprocessSuffix) "matrix.c"

$(IntermediateDirectory)/newtonmin$(ObjectSuffix): newtonmin.c $(IntermediateDirectory)/newtonmin$(DependSuffix)
	$(CC) $(SourceSwitch) "C:/Documents and Settings/USER/My Documents/Downloads/codelite/ctest/newtonmin/newtonmin.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/newtonmin$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/newtonmin$(DependSuffix): newtonmin.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/newtonmin$(ObjectSuffix) -MF$(IntermediateDirectory)/newtonmin$(DependSuffix) -MM "newtonmin.c"

$(IntermediateDirectory)/newtonmin$(PreprocessSuffix): newtonmin.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/newtonmin$(PreprocessSuffix) "newtonmin.c"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) $(IntermediateDirectory)/main$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/main$(DependSuffix)
	$(RM) $(IntermediateDirectory)/main$(PreprocessSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(DependSuffix)
	$(RM) $(IntermediateDirectory)/matrix$(PreprocessSuffix)
	$(RM) $(IntermediateDirectory)/newtonmin$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/newtonmin$(DependSuffix)
	$(RM) $(IntermediateDirectory)/newtonmin$(PreprocessSuffix)
	$(RM) $(OutputFile)
	$(RM) $(OutputFile).exe
	$(RM) "../.build-debug/newtonmin"


