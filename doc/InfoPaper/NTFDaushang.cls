%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}
%
\newcommand{\NTFD@name}{NTFDaushang}
\newcommand{\NTFD@date}{2011/05/17}
\newcommand{\NTFD@version}{0.2}
\newcommand{\NTFD@msg}{Aushang im Design des NTFD}
\newcommand{\NTFD@ack}{In recognition of the work of Thomas Benkert}
\newcommand{\NTFD@sep}{############################################################}
%
\typeout{\NTFD@sep}
\ProvidesClass{\NTFD@name}[\NTFD@date\space v.\NTFD@version\space Felix Dietzsch]
%
\typeout{\space\space\space\space\space\space\space\space\space\space\space\space\space\space\space\space\NTFD@msg}
\typeout{\space\space\space\space\space\space\space\space\space\space\space\space\space\space\space\space\NTFD@ack}
\typeout{\NTFD@sep}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grundklasse laden
\LoadClass[a4paper]{scrartcl}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pakete für die Optionen
\RequirePackage{xkeyval}
\RequirePackage{ifthen}
\RequirePackage{xifthen}
\RequirePackage{etoolbox}
% weitere Pakete
\RequirePackage[]{scrpage2}
%\RequirePackage{lastpage}
\RequirePackage[%
	paper=a4paper,
	inner=78.81409pt,%\the\NTFDbisT,% wg. druck ohne bausteine
	outer=21.9mm,%\the\NTFDabMitteLogo,% wg. druck ohne bausteine
	top=18.3mm,
	bottom=30mm,%17.3mm
	headsep=4.3ex,
	footskip=20mm,
	includehead
	]{geometry}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neue Befehle für Optionen
\newcommand{\NTFDBulletin@parskip}{}
\newcommand{\NTFDBulletin@fontsize}{}
\newcommand{\NTFDBulletin@sprache}{}
\newboolean{NTFDBulletin@futura}\setboolean{NTFDBulletin@futura}{false}%
\newboolean{NTFDBulletin@german}\setboolean{NTFDBulletin@german}{true}%
\newboolean{NTFDBulletin@druck}\setboolean{NTFDBulletin@druck}{false}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optionen deklarieren
%
% für übergabe an Koma
\DeclareOptionX{parskip}{\renewcommand\NTFDBulletin@parskip{#1}}
%
\DeclareOptionX{fontsize}{\renewcommand\NTFDBulletin@fontsize{#1}}
%
\DeclareOptionX{futura}[true]{%
	\XKV@cc*+[\val\nr]{#1}{true,false}{%
	\ifcase\nr\relax
		\setboolean{NTFDBulletin@futura}{true}%
	   \or
		\setboolean{NTFDBulletin@futura}{false}%
	\fi
	}{\PackageWarning{\NTFDBulletin@name}{the input 'futura=\val' was incorrect and was ignored.}}%
}
%
\DeclareOptionX{stil}[epost]{%
	\XKV@cc*+[\val\nr]{#1}{epost,druck}{%
	\ifcase\nr\relax
		\setboolean{NTFDBulletin@druck}{false}%
	   \or
		\setboolean{NTFDBulletin@druck}{true}%
	\fi
	}{\PackageWarning{\NTFDBulletin@name}{the input 'stil=\val' was incorrect and was ignored.}}%
}
\DeclareOptionX{sprache}[ngerman]{%
	\XKV@cc*+[\val\nr]{#1}{ngerman,english}{%
	\ifcase\nr\relax
		\renewcommand{\NTFDBulletin@sprache}{ngerman}%
	   \or
		\renewcommand{\NTFDBulletin@sprache}{english}%
		\setboolean{NTFDBulletin@german}{false}%
	\fi
	}{\ClassWarning{\NTFDBulletin@name}{the input 'sprache=\val' was incorrect and was ignored.}}%
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
\DeclareOptionX*{\ClassWarning{\NTFDBulletin@name}{Unknown option ignored: \CurrentOption}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optionen ausführen
% set options
\ExecuteOptionsX{%
	fontsize=11pt,
	parskip=half,
	futura=false,
	sprache=english,
	stil=epost
}
% alle Optionen abarbeiten
\ProcessOptionsX\relax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testen auf Futura
\ifthenelse{\boolean{NTFDBulletin@futura}}{%
	\RequirePackage{TUBAFfutura}%
   }{%
}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testen auf stil
%\ifthenelse{\boolean{NTFDBulletin@druck}}{%
%   }{%
% is mandatory for the inclusion of the header and footer
\RequirePackage{NTFDbausteine}
%}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nachtrügliches Setzen der Klassenoptionen
\KOMAoptions{%
	fontsize=\NTFDBulletin@fontsize,
	parskip=\NTFDBulletin@parskip
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sprachanpassungen
\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
	\RequirePackage[\NTFDBulletin@sprache]{babel}%
   }{%
	\RequirePackage[ngerman,\NTFDBulletin@sprache]{babel}%
}%
\input{NTFDaushang.ldf}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neue Befehle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name der Fakultaet
\def\NTFDFakultaet#1{\def\insertNTFDFakultaet{#1}}
	\def\insertNTFDFakultaet{}
\NTFDFakultaet{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name des Institutes
\def\NTFDInstitut#1{\def\insertNTFDInstitut{#1}}
	\def\insertNTFDInstitut{}
\NTFDInstitut{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name des Lehrstuhls
\def\NTFDLehrstuhl#1{\def\insertNTFDLehrstuhl{#1}}
	\def\insertNTFDLehrstuhl{}
\NTFDLehrstuhl{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Studenten
\def\NTFDStudents{\NTFDBulletin@Students}
\long\def\NTFDBulletin@Students[#1][#2][#3]#4{
    \def\NTFDBulletin@temp{#4}
	\ifdefempty{\NTFDBulletin@temp}{
		\def\NTFDStudentName{}
		\def\insertNTFDStudentOne{}
		\def\insertNTFDStudentTwo{}
		\def\insertNTFDStudentThree{}
		\def\insertNTFDStudentFour{}
		}{
		\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
			\def\NTFDStudentName{Studenten:}
			}{
			\def\NTFDStudentName{Students:}
			}
        \def\insertNTFDStudentOne{#4}%
		\def\insertNTFDStudentTwo{}
		\def\insertNTFDStudentThree{}
		\def\insertNTFDStudentFour{}
		}
    \def\NTFDBulletin@temp{#3}
	\ifdefempty{\NTFDBulletin@temp}{}
		{\def\insertNTFDStudentTwo{#3}}
    \def\NTFDBulletin@temp{#2}
	\ifdefempty{\NTFDBulletin@temp}{}
		{\def\insertNTFDStudentThree{#2}}
    \def\NTFDBulletin@temp{#1}
	\ifdefempty{\NTFDBulletin@temp}{}
		{\def\insertNTFDStudentFour{#1}}
}
\NTFDStudents[][][]{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gutachter
\def\NTFDGutachter{\@dblarg\NTFDBulletin@Gutachter}
\long\def\NTFDBulletin@Gutachter[#1]#2{%
    \def\NTFDBulletin@temp{#1}%
    \ifx\NTFDBulletin@temp\@empty%
		\def\NTFDGutachtername{}
        \def\insertNTFDGutachter{}%
    \else%
		\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
			\def\NTFDGutachtername{Gutachter:}
			}{
			\def\NTFDGutachtername{Assessor:}
			}
        \def\insertNTFDGutachter{#1}%
    \fi%
}%
\NTFDGutachter{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Betreuer
\def\NTFDBetreuer{\NTFDBulletin@Betreuer}
\long\def\NTFDBulletin@Betreuer[#1][#2]#3{%
    \def\NTFDBulletin@temp{#3}%
    \ifx\NTFDBulletin@temp\@empty%
		\def\NTFDBetreuername{}
        \def\insertNTFDBetreuerE{}%
    \else%
		\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
			\def\NTFDBetreuername{Betreuer:}
			}{
			\def\NTFDBetreuername{Advisor:}
		}
        \def\insertNTFDBetreuerE{#3}
    \fi%
    \def\NTFDBulletin@temp{#2}%
    \ifx\NTFDBulletin@temp\@empty%
        \def\insertNTFDBetreuerZ{}%
    \else%
		\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
			\def\NTFDBetreuername{Betreuer:}
			}{
			\def\NTFDBetreuername{Advisors:}
		}
        \def\insertNTFDBetreuerZ{#2}
    \fi%
    \def\NTFDBulletin@temp{#1}%
    \ifx\NTFDBulletin@temp\@empty%
        \def\insertNTFDBetreuerT{}%
    \else%
        \def\insertNTFDBetreuerT{#1}
    \fi%
}%
\NTFDBetreuer[][]{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name des Bulletin
\def\NTFDBulletinName#1{\def\insertNTFDBulletinName{#1}}
	\def\insertNTFDBulletinName{}
\NTFDBulletinName{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Titel der Arbeit
\def\NTFDTitel{\@dblarg\NTFDBulletin@Titel}
\long\def\NTFDBulletin@Titel[#1]#2{%
  \def\NTFDBulletin@temp{#2}%
  \ifx\NTFDBulletin@temp\@empty
    \def\insertNTFDTitel{}
  \else
    \def\insertNTFDTitel{#2}%
  \fi
 \def\insertshortNTFDTitel{#1}
}
\NTFDTitel{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor
\def\NTFDAutor{\@dblarg\NTFDBulletin@Autor}
\long\def\NTFDBulletin@Autor[#1]#2{%
  \def\NTFDBulletin@temp{#2}%
  \ifx\NTFDBulletin@temp\@empty
    \def\insertNTFDAutor{}
  \else
    \def\insertNTFDAutor{#2}%
  \fi
 \def\insertshortNTFDAutor{#1}}
\NTFDAutor{}
\newcommand{\NTFDAutorprefix}{vorgelegt~von}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Datum
\def\NTFDDatum{\@dblarg\NTFDBulletin@Datum}
\long\def\NTFDBulletin@Datum[#1]#2{%
  \def\NTFDBulletin@temp{#2}%
  \ifx\NTFDBulletin@temp\@empty
    \def\insertNTFDDatum{}
  \else
    \def\insertNTFDDatum{#2}%
  \fi
 \def\insertshortNTFDDatum{#1}}
\NTFDDatum{\today}
\newcommand{\NTFDDatumprefix}{vom}
\newcommand{\shortNTFDDatumprefix}{(}
\newcommand{\shortNTFDDatumsuffix}{)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F"ullw"orter f"ur die Titelseite
%\newcommand{\NTFDBulletinName}{}
%
\def\NTFDBulletinSuffixE#1{\def\insertNTFDBulletinSuffixE{#1}}
	\def\insertNTFDBulletinSuffixE{}
\NTFDBulletinSuffixE{}
\def\NTFDBulletinSuffixZ#1{\def\insertNTFDBulletinSuffixZ{#1}}
	\def\insertNTFDBulletinSuffixZ{}
\NTFDBulletinSuffixZ{}
\newcommand{\NTFDseitenname}{Seite}
\newcommand{\NTFDseitentrenner}{~von~}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Umgebung Bulletin
%\newenvironment{NTFDBulletin}{%
%	\maketitle%
%	\noindent\ignorespaces%
%   }{%
%	\par\vspace*{4\baselineskip}
%	\begin{tabular}[tbh]{ll}
%	   \noindent{\bfseries\NTFDBetreuername}&\insertNTFDBetreuerE\\
%	   \noindent&\insertNTFDBetreuerZ\\
%	   \noindent&\insertNTFDBetreuerT\\
%       \noindent{\bfseries\NTFDGutachtername}&\insertNTFDGutachter
%	\end{tabular}
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neue Befehle Ende 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% und hier die Redefinition der Titelseite
%
\newcommand{\NTFD@yoffset}{0}
\newcommand{\NTFD@xoffset}{0}
\newcommand{\NTFDbausteine@breite}{}
\newcommand{\NTFDbausteine@basislaenge}{}
\newcommand{\NTFDbausteine@basiseinheit}{mm}
\newcommand{\NTFDbausteine@tmpa}{}
\newcommand{\NTFDbausteine@scale}{1}
\newlength\NTFD@temp

%\TPshowboxestrue % show boxes of textblocks
 
\renewcommand*{\titlepagestyle}{scrheadings}
\renewcommand*{\maketitle}{%
	\thispagestyle{scrheadings}
	\clearscrheadfoot
	\ifthenelse{\boolean{NTFDBulletin@druck}}{%
       		\thispagestyle{empty}
		\vspace*{0pt}%
	 	}{%
       		\thispagestyle{empty}
		%\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{%
			\ifthenelse{\lengthtest{\parskip=\z@}}{%
               			 \begin{textblock*}{\paperwidth}(-11pt,0pt)% 11pt sonst overfull hbox
           		}{% 
                			\begin{textblock*}{\paperwidth}(0pt,0pt)
        			}   
                			\ifthenelse{\boolean{NTFDBulletin@druck}}{%
                   		}{% 
                        			\TUBAFkopf%
                			}
				\end{textblock*}
				%##############################
				% NTFD Logo
				\setlength{\NTFD@temp}{37.337mm}% Korrektur, sonst ovefull hbox
				\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
				
				\setlength{\TUBAFbausteine@tmp}{\TUBAFbausteine@bausteinelaenge}
				\FPset\TUBAFbausteine@basislaenge{\number\TUBAFbausteine@tmp}
				
				\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge				
				\FPround\NTFDbausteine@breite\NTFDbausteine@tmpa5
				
				\FPset{\myX}{5.92}
				\FPset{\myY}{4.75}
				\FPset{\NTFDbausteine@scale}{1.06}
				
				\FPsub\NTFD@xoffset\TUBAFbausteine@textboxbreite\TUBAFbausteineRandRiv
				\FPsub\NTFD@xoffset\NTFD@xoffset\NTFDbausteine@breite
				\FPsub\NTFD@xoffset\NTFD@xoffset\myX
				
				\FPset{\NTFD@yoffset}{\TUBAFbausteineRandOiii}
				\FPadd\NTFD@yoffset\NTFD@yoffset\myY
				
				\begin{textblock}{\NTFDbausteine@breite}(\NTFD@xoffset,\NTFD@yoffset)
					\renewcommand{\NTFDlogofgLight}{NTFDrot}
					\renewcommand{\NTFDlogofgDark}{NTFDgrau}
					\begin{tikzpicture}[scale=\NTFDbausteine@scale]
						\NTFDLogoeinzeln
					\end{tikzpicture}
				\end{textblock}
				% Schriftzug Professur
				\setlength{\NTFD@temp}{99.206mm}% Korrektur, sonst ovefull hbox
				\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
				
				\setlength{\TUBAFbausteine@tmp}{\TUBAFbausteine@bausteinelaenge}
				\FPset\TUBAFbausteine@basislaenge{\number\TUBAFbausteine@tmp}
				
				\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge				
				\FPround\NTFDbausteine@breite\NTFDbausteine@tmpa5
				
				% Y Offset
				\FPset{\NTFD@yoffset}{\TUBAFbausteineRandOiv}
				\FPadd\NTFD@yoffset\NTFD@yoffset\myY
				% X Offset
				\FPset{\myX}{0.6}
				\FPset{\NTFD@xoffset}{\TUBAFbausteineRandLi}
				\FPsub\NTFD@xoffset\NTFD@xoffset\myX
				% Scaling
				\FPset{\NTFDbausteine@scale}{1.0}
				
				\begin{textblock}{\NTFDbausteine@breite}(\NTFD@xoffset,\NTFD@yoffset)
					\renewcommand{\NTFDlogofgLight}{NTFDrot}
					\renewcommand{\NTFDlogofgDark}{NTFDgrau}
					\begin{tikzpicture}[scale=\NTFDbausteine@scale]
						\ifthenelse{\equal{\NTFDBulletin@sprache}{ngerman}}{
							\NTFDAushangKopfSchriftD
						}{
							\NTFDAushangKopfSchriftE
						}
					\end{tikzpicture}
				\end{textblock}				
			%}
		}
		
		\setlength{\NTFD@temp}{1\textwidth}% Korrektur, sonst ovefull hbox
		\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
		\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge
		\FPround\NTFDbausteine@breite\NTFDbausteine@tmpa5
		% Y Offset	
		\FPset{\myY}{10}
		\FPset{\NTFD@yoffset}{\TUBAFbausteineRandOii}
		\FPadd\NTFD@yoffset\NTFD@yoffset\myY
		% X position
		\FPset{\myX}{5.5}
		\FPset{\NTFD@xoffset}{\TUBAFbausteineRandRi}
		\FPadd\NTFD@xoffset\NTFD@xoffset\myX

		\begin{textblock}{\NTFDbausteine@breite}(\NTFD@xoffset,\NTFD@yoffset)
			\sffamily%
			{\centering\bfseries\LARGE\insertNTFDTitel\vskip1.8ex}
			%\ifdefempty{\insertNTFDBulletinName}{}{\centering\bfseries\LARGE\insertNTFDBulletinName\vskip1.5ex}%
			%\ifdefempty{\insertNTFDBulletinSuffixE}{}{\normalsize{\insertNTFDBulletinSuffixE}\vskip1.5ex}
			%{\normalsize{\insertNTFDBulletinSuffixZ}\vskip2ex}
		\end{textblock}
		% Fusszeile
		% Breite der Fuszeile
		\setlength{\NTFD@temp}{182.483mm}% Korrektur, sonst ovefull hbox
		\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
		\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge
		\FPround\NTFDbausteine@breite\NTFDbausteine@tmpa5
		
		% Y position
		\setlength{\NTFD@temp}{\paperheight}% Korrektur, sonst ovefull hbox
		\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
		\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge
		\FPround\NTFDbausteine@tmpa\NTFDbausteine@tmpa5
		\FPset{\myY}{15.4}
		\FPsub\NTFDbausteine@tmpa\NTFDbausteine@tmpa\myY
		
		\ifthenelse{\boolean{NTFDBulletin@druck}}{
		}{%
			\begin{textblock}{\NTFDbausteine@breite}(\NTFD@xoffset,\NTFDbausteine@tmpa)
				\renewcommand{\NTFDFootFG}{NTFDgrau}
				\begin{tikzpicture}
					\NTFDFusszeile
				\end{tikzpicture}
			\end{textblock}
		}
% Rest des Dokumentes	
\FPset{\myY}{0}
\FPadd\NTFD@yoffset\NTFD@yoffset\myY
\vspace*{\NTFD@yoffset\NTFDbausteine@basiseinheit}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gestaltung des restlichen Dokumentes
\AtBeginDocument{
	\sffamily%
	\pagestyle{bulletin}
}
% Define own page style
\newcommand{\ps@bulletin}{
	\renewcommand{\@oddhead}{}% header blank
	\renewcommand{\@evenhead}{}% header blank
	\renewcommand{\@oddfoot}{
				\renewcommand{\NTFDFootFG}{NTFDgrau}
				\begin{tikzpicture}
				\put(-5,-28.3){
					\NTFDFusszeile
					}
				\end{tikzpicture}
	}%
	\renewcommand{\@evenfoot}{
				\renewcommand{\NTFDFootFG}{NTFDgrau}
				\begin{tikzpicture}
				\put(-100,-28.5){
					\NTFDFusszeile
					}
				\end{tikzpicture}
	}%
}
%\renewcommand{\ps@bulletin}{
%	\renewcommand{\@oddhead}{%     header blank
%	\renewcommand{\@evenhead}{%    header blank
%	\renewcommand{\@oddfoot}{}%     
%    \renewcommand{\@evenfoot}{}%     
%	  }
%	}}
\newenvironment{NTFDBulletin}{
	\maketitle
	\noindent\ignorespaces%
}{
	\par\vspace*{4\baselineskip}
	\begin{tabular}[tbh]{ll}
		{\bfseries\NTFDStudentName} & \insertNTFDStudentOne\\
		\ifdefempty{\insertNTFDStudentTwo}{}{\noindent & \insertNTFDStudentTwo\\}
		\ifdefempty{\insertNTFDStudentThree}{}{\noindent & \insertNTFDStudentThree\\}
		\ifdefempty{\insertNTFDStudentFour}{}{\noindent & \insertNTFDStudentFour\\}
	\noindent{\bfseries\NTFDBetreuername}&\insertNTFDBetreuerE\\
	\ifdefempty{\insertNTFDBetreuerZ}{}{\noindent&\insertNTFDBetreuerZ\\}
	\ifdefempty{\insertNTFDBetreuerT}{}{\noindent&\insertNTFDBetreuerT\\}
	\noindent{\bfseries\NTFDGutachtername}&\insertNTFDGutachter
	\end{tabular}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schriftanpassungen
\setkomafont{captionlabel}{\footnotesize\sffamily\bfseries}
\setkomafont{caption}{\footnotesize}
\setcapindent{0pt}
\setkomafont{pageheadfoot}{\normalfont\footnotesize\sffamily}
\setkomafont{pagenumber}{\normalfont\normalsize\sffamily}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anstatt des frame-Befehls; arbeitet in allen Modi
\long\def\NTFDumrandung#1{%
  \leavevmode
  \hbox{%
    \hskip-\@wholewidth
    \vbox{%
      \vskip-\@wholewidth
      \hrule \@height\@wholewidth
      \hbox{%
        \vrule\@width\@wholewidth
        #1%
        \vrule\@width\@wholewidth}%
      \hrule\@height\@wholewidth
      \vskip-\@wholewidth}%
    \hskip-\@wholewidth}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\endinput
