% copyright 2008-2009 by Thomas Benkert/U.Haack
% copyright 2009-2011 by Thomas Benkert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dieses Werk ist unter einem Creative Commons 
% Namensnennung - Keine Bearbeitung 3.0 Deutschland Lizenzvertrag lizenziert. 
% Um die Lizenz anzusehen, gehen Sie bitte zu http://creativecommons.org/licenses/by-nd/3.0/de/ 
% oder schicken Sie einen Brief an 
% Creative Commons, 171 Second Street, Suite 300, San Francisco, California 94105, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}
%
\newcommand{\NTFDbrief@name}{NTFDbrief}
\newcommand{\NTFDbrief@date}{2011/11/09}
\newcommand{\NTFDbrief@version}{1.5.0}
%
\ProvidesClass{\NTFDbrief@name}[\NTFDbrief@date\space v.\NTFDbrief@version\space Thomas Benkert]
%
\typeout{\@spaces\@spaces\@spaces\@spaces eine LaTeX Briefklasse im Corporate Design^^J
	\@spaces\@spaces\@spaces\space\space\space der TU Bergakademie Freiberg}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scrlttr2 laden
\LoadClass[%
a4paper,%
pagenumber=off,%
refline=narrow]{scrlttr2}[2009/04/03]%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pakete einbinden
\RequirePackage{xkeyval}%
\RequirePackage{ifthen}%
\RequirePackage{textcomp}%
\RequirePackage{lastpage}%
\RequirePackage{scrpage2}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neue Befehle
\newcommand{\NTFDbrief@parskip}{}
\newcommand{\NTFDbrief@fontsize}{}
\newcommand{\NTFDbrief@hausschrift}{}
\newboolean{NTFDbrief@hausschrift}\setboolean{NTFDbrief@hausschrift}{true}%
\newcommand{\NTFDbrief@stil}{normal}
\newcommand{\NTFDbrief@RetouradresseHkor}{0cm}
\newcommand{\NTFDbrief@sprache}{}
\newcommand{\NTFDbrief@funktion}{}
\newcommand{\NTFDbrief@zusatztext}{normal}
%
\newcommand{\NTFDbrief@bausteinedvipdfm}{false}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Klassenoptionen
%
\DeclareOptionX{parskip}{\renewcommand\NTFDbrief@parskip{#1}}
%
\DeclareOptionX{fontsize}[11pt]{%
	\XKV@cc*+[\val\nr]{#1}{10pt,11pt,12pt}{%
	\ifcase\nr\relax
		\renewcommand{\NTFDbrief@fontsize}{10pt}%
	   \or
		\renewcommand{\NTFDbrief@fontsize}{11pt}%
	   \or
		\renewcommand{\NTFDbrief@fontsize}{12pt}%
	\fi    
	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'fontsize=\val' is incorrect and ignored}}%
}
%
\DeclareOptionX{futura}[true]{%
	\XKV@cc*+[\val\nr]{#1}{true,false}{%
	\ifcase\nr\relax
		\setboolean{NTFDbrief@hausschrift}{true}%
		\ClassWarningNoLine{\NTFDbrief@name}{The option 'futura' is obsolete.^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\space\space\space
			Please use 'hausschrift' instead^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces}%
	   \or
		\setboolean{NTFDbrief@hausschrift}{false}%
		\ClassWarningNoLine{\NTFDbrief@name}{The option 'futura' is obsolete.^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\space\space\space
			Please use 'hausschrift' instead^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces}%
	\fi    
	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'futura=\val' was incorrect and was ignored.^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\space\space\space
			The option 'futura' is obsolete.^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\space\space\space
			Please use 'hausschrift' instead^^J
			\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces\@spaces}%
	}%
}
%
\DeclareOptionX{hausschrift}[true]{%
	\ifthenelse{\equal{#1}{false}}{%
		\setboolean{NTFDbrief@hausschrift}{false}%
	   }{%
	   	\renewcommand\NTFDbrief@hausschrift{#1}%
		\setboolean{NTFDbrief@hausschrift}{true}%
	}{}%\ClassWarningNoLine{\NTFDarbeiten@name}{The input 'hausschrift=\val' is incorrect and ignored}}%
}
%
\DeclareOptionX{stil}[normal]{%
	\XKV@cc*+[\val\nr]{#1}{normal,intern,epostsimpel,epost,wagner,ernst,druck}{%
	\ifcase\nr\relax
		\renewcommand{\NTFDbrief@stil}{normal}%
	   \or
		\renewcommand{\NTFDbrief@stil}{intern}%
	   \or
		\renewcommand{\NTFDbrief@stil}{epostsimpel}%
	   \or
		\renewcommand{\NTFDbrief@stil}{epost}%
	   \or
		\renewcommand{\NTFDbrief@stil}{wagner}%
	   \or
		\renewcommand{\NTFDbrief@stil}{ernst}%
	   \or
		\renewcommand{\NTFDbrief@stil}{druck}%
	\fi    
	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'stil=\val' is incorrect and ignored}}%
}
%
\DeclareOptionX{rahkor}{\renewcommand{\NTFDbrief@RetouradresseHkor}{1.5cm}}
%
\DeclareOptionX{sprache}[ngerman]{%
	\XKV@cc*+[\val\nr]{#1}{nobabel,ngerman,german,american,english,french,norsk,polish}{%
	\ifcase\nr\relax
		\renewcommand{\NTFDbrief@sprache}{nobabel}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{ngerman}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{german}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{american}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{english}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{french}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{norsk}%
	   \or
		\renewcommand{\NTFDbrief@sprache}{polish}%
	\fi
	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'sprache=\val' is incorrect and ignored}}%
}
%
\DeclareOptionX{funktion}[beides]{%
	\XKV@cc*+[\val\nr]{#1}{oben,unten,beides}{%
	\ifcase\nr\relax
		\renewcommand{\NTFDbrief@funktion}{oben}%
	   \or
		\renewcommand{\NTFDbrief@funktion}{unten}%
	   \or
		\renewcommand{\NTFDbrief@funktion}{beides}%
	\fi    
	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'funktion=\val' is incorrect and ignored}}%
}
%
\DeclareOptionX{zusatz}{\renewcommand{\NTFDbrief@zusatztext}{#1}}
%%
%\DeclareOptionX{zusatz}[normal]{%
%	\XKV@cc*+[\val\nr]{#1}{hell,normal,dunkel,dunkler,schwarz}{%
%	\ifcase\nr\relax
%		\renewcommand{\NTFDbrief@zusatztext}{NTFDzusatzhell}%
%	   \or
%		\renewcommand{\NTFDbrief@zusatztext}{NTFDzusatz}%
%	   \or
%		\renewcommand{\NTFDbrief@zusatztext}{NTFDzusatzdunkel}%
%	   \or
%		\renewcommand{\NTFDbrief@zusatztext}{NTFDzusatzdunkler}%
%	   \or
%		\renewcommand{\NTFDbrief@zusatztext}{NTFDschwarz}%
%	\fi
%	}{\ClassWarningNoLine{\NTFDbrief@name}{the input 'zusatz=#1' was incorrect and was ignored.}}%
%}
%
\DeclareOptionX{dvipdfm}[true]{\renewcommand\NTFDbrief@bausteinedvipdfm{#1}}
%
\DeclareOptionX*{\ClassWarningNoLine{\NTFDbrief@name}{Unknown option ignored: \CurrentOption}}
%
\ExecuteOptionsX{%
	hausschrift=false,
	stil=normal,
	fontsize=11pt,
	parskip=half,
	sprache=ngerman,
	funktion=beides,
}
\ProcessOptionsX\relax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nachträgliches Setzen der Klassenoptionen
\KOMAoptions{%
	fontsize=\NTFDbrief@fontsize,%
	parskip=\NTFDbrief@parskip%
	}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testen auf Hausschrift
\ifthenelse{\boolean{NTFDbrief@hausschrift}}{%
	\RequirePackage[staerke=\NTFDbrief@hausschrift]{NTFDhausschrift}
	}{%
	\renewcommand{\familydefault}{cmss}%
}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testen auf stil und Laden der Pakete
%\ifthenelse{\equal{\NTFDbrief@stil}{intern}\OR\equal{\NTFDbrief@stil}{epostsimpel}\OR\equal{\NTFDbrief@stil}{epost}}{%
	%\RequirePackage[dvipdfm=\NTFDbrief@bausteinedvipdfm]{TUBAFbausteine}%
   %}{%
%}
%\RequirePackage{NTFDbausteinelaengen}%
\RequirePackage{NTFDbausteine}
\RequirePackage[%
	left=\TUBAFbausteineRandLi\TUBAFbausteine@bausteineeinheit,
	right=\TUBAFbausteineRandRi\TUBAFbausteine@bausteineeinheit,
	top=19.4mm,
	bottom=14.2pt,
	includeheadfoot
	]{geometry}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zusatz setzen
\ifthenelse{\equal{\NTFDbrief@stil}{intern}}{%\OR\equal{\NTFDbrief@stil}{epost}}{%
	\ifthenelse{\equal{\NTFDbrief@zusatztext}{hell}}{%
		\renewcommand{\NTFDzusatztext}{NTFDzusatzhell}}{%
		\ifthenelse{\equal{\NTFDbrief@zusatztext}{normal}}{%
			\renewcommand{\NTFDzusatztext}{NTFDzusatz}}{%
			\ifthenelse{\equal{\NTFDbrief@zusatztext}{dunkel}}{%
				\renewcommand{\NTFDzusatztext}{NTFDzusatzdunkel}}{%
				\ifthenelse{\equal{\NTFDbrief@zusatztext}{dunkler}}{%
					\renewcommand{\NTFDzusatztext}{NTFDzusatzdunkler}}{%
					\ifthenelse{\equal{\NTFDbrief@zusatztext}{schwarz}}{%
						\renewcommand{\NTFDzusatztext}{NTFDschwarz}}{%
				\ClassWarningNoLine{\NTFDbrief@name}{the input 'zusatz=\NTFDbrief@zusatztext' was incorrect and was ignored.}%
	}}}}}%
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sprachanpassungen
\ifthenelse{\equal{\NTFDbrief@sprache}{nobabel}}{%
	\ClassWarningNoLine{\NTFDbrief@name}{Package babel not loaded.\MessageBreak
							Fall back to default language.}
   }{%
	\ifthenelse{\equal{\NTFDbrief@sprache}{ngerman}}{%
		\RequirePackage[\NTFDbrief@sprache]{babel}%
	   }{%
		\RequirePackage[ngerman,\NTFDbrief@sprache]{babel}%
	}%
	\input{NTFDbrief.ldf}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudolaengen anpassen
\@setplength{toaddrvpos}{45mm}%
\@setplength{toaddrhpos}{\TUBAFbausteineRandLi\TUBAFbausteine@bausteineeinheit}%
\@setplength{refhpos}{\TUBAFbausteineRandLi\TUBAFbausteine@bausteineeinheit}%
\@setplength{backaddrheight}{5mm}%
\@setplength{toaddrheight}{45mm}%
\@setplength{toaddrwidth}{85mm}%
\@setplength{toaddrindent}{0mm}%
\@setplength{specialmailrightindent}{0mm}% rechte Einrückung Versandart
\ifthenelse{\equal{\NTFDbrief@RetouradresseHkor}{0cm}}{}{%
	\@addtoplength{toaddrhpos}{-\NTFDbrief@RetouradresseHkor}%
	\@addtoplength{toaddrindent}{\NTFDbrief@RetouradresseHkor}%
}
\@setplength{tfoldmarklength}{3mm}%	obere Falsmarke Länge
\@setplength{tfoldmarkvpos}{105mm}%	obere Faltmarke Abstand von oberer Papierkante
\@setplength{pfoldmarklength}{1.6\useplength{tfoldmarklength}}% Lochermarke
\@setplength{bfoldmarkvpos}{210mm}%	untere Faltmarke Länge
\@setplength{bfoldmarklength}{\useplength{tfoldmarklength}}% untere Faltmarke 
\@setplength{foldmarkhpos}{3mm}%	horizonaler Abstand Marken
\@setplength{lochpos}{0mm}%	location
\@setplength{locvpos}{103mm}%	location-Feld vertikaler Abstand
\@setplength{sigbeforevskip}{3\baselineskip}% Abstand Grußformel/Signatur
\@setplength{firstheadvpos}{42.2pt}%
\@setplength{firstheadwidth}{140mm}%gemessen
\@addtoplength{subjectaftervskip}{-6pt}
\@setplength{firstfootvpos}{800.9pt}%{815.8pt}%
\@setplength{firstfoothpos}{74.5pt}%{815.8pt}%
\ifthenelse{\equal{\NTFDbrief@fontsize}{10pt}}{%
	\@setplength{firstfootvpos}{816.2pt}}{}
\@setplength{firstfootwidth}{\paperwidth}%
\setlength{\footskip}{22.47482pt}%1.75ex+5mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stellen für Zeilenumbrüche in URLs
\def\UrlBreaks{\do\.\do\@\do\\\do\/\do\!\do\_\do\|\do\;\do\>\do\]\do\)\do\,\do\?\do\'\do\+\do\=\do\#\do\-}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seitenzahlen für mehrere Briefe ermöglichen, gibt Warnung beim ersten Kompilieren
\newcounter{letter}%
\AtBeginLetter{\stepcounter{letter}}%
\renewcommand*{\stopletter}{\par\label{lastpage@\theletter}}%
\newcommand*{\pagesatthisletter}{\pageref{lastpage@\theletter}}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die offsets einstellen in mm
\newcommand{\NTFDbriefxoffset}{35.01}% gemessen
\newcommand{\NTFDbriefyoffset}{-11.37}% gemessen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Briefkopf erste Seite
\newcommand{\NTFD@yoffset}{0}
\newcommand{\NTFD@xoffset}{0}
\newcommand{\NTFDbausteine@breite}{}
\newcommand{\NTFDbausteine@basislaenge}{}
\newcommand{\NTFDbausteine@basiseinheit}{mm}
\newcommand{\NTFDbausteine@tmpa}{}
\newcommand{\NTFDbausteine@scale}{1}
\newlength\NTFD@temp


\setkomavar{firsthead}{%
	\ifthenelse{\equal{\NTFDbrief@stil}{epost}}{%
		\ifthenelse{\lengthtest{\parskip=\z@}}{%
               			 \begin{textblock*}{\paperwidth}(-11pt,0pt)% 11pt sonst overfull hbox
           		}{% 
                			\begin{textblock*}{\paperwidth}(0pt,0pt)
        			}   
						\ifthenelse{\equal{\NTFDbrief@stil}{druck}}{%
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
		   }{}%
		%\ifthenelse{\equal{\NTFDbrief@stil}{epostsimpel}\OR\equal{\NTFDbrief@stil}{epost}}{%
			%\TUBAFrahmen[fg=weiss]%
		   %}{%
		%}%
		%\ifthenelse{\lengthtest{\parskip=\z@}}{%
			%\ifthenelse{\equal{\NTFDbrief@stil}{ernst}}{%
				%\hskip-22pt}{\hskip-33pt}%
		   %}{%
			%\hskip-22pt%
		%}% \hskip gemessen
		%\begin{minipage}[t]{\useplength{firstheadwidth}}
%%			\vskip11.8mm% ebenso
			%\vskip14mm% korregiert wg. fehlerhafter Briefbögen und Druckereinzugs
			%\ifkomavarempty{fakultaet}{}{\usekomavar*{fakultaet}\usekomavar{fakultaet}}%
			%\ifkomavarempty{institut}{}{\newline\usekomavar*{institut}\usekomavar{institut}}%
			%\ifkomavarempty{lehrstuhl}{}{\newline\usekomavar*{lehrstuhl}\usekomavar{lehrstuhl}}%
			%\ifkomavarempty{funktion}{}{%
				%\ifthenelse{\equal{\NTFDbrief@funktion}{beides}\OR%
					      %\equal{\NTFDbrief@funktion}{oben}}{%
					%\newline\usekomavar*{funktion}\usekomavar{funktion}}{}%
			%}%
		%\end{minipage}%
		%\vskip-5mm% Warnung (wenn alles ausgefüllt) unterdrücken
	}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mit location Informationsblock nach DIN 676:1995-05 Bild 5 definieren
\newcommand*{\NTFDrefnames}{}
\newcommand*{\NTFDrefs}{}
\AtBeginLetter{\setkomavar{location}{%
	\ifkomavarempty{yourref}{}{%
		\l@addto@macro\NTFDrefnames{\usekomavar*{yourref} &}}%: & \usekomavar{yourref} \\}}%
	\ifkomavarempty{yourmail}{}{%
		\l@addto@macro\NTFDrefnames{\usekomavar*{yourmail} &}}%: & \usekomavar{yourmail} \\}}%
	\ifkomavarempty{myref}{}{%
		\l@addto@macro\NTFDrefnames{\usekomavar*{myref} &}}%: & \usekomavar{myref} \\}}%
	\ifkomavarempty{unserschreiben}{}{%
		\l@addto@macro\NTFDrefnames{\usekomavar*{unserschreiben}}}% & \usekomavar{unserschreiben} \\}}%
	\l@addto@macro\NTFDrefnames{\\}%
	%\ifkomavarempty{fromname}{}{%
		%\l@addto@macro\NTFDlocation{\usekomavar*{fromname}: & \usekomavar{fromname} \\}}%
	%\ifkomavarempty{durchwahl}{}{%
		%\l@addto@macro\NTFDlocation{\usekomavar*{durchwahl} & \NTFDTelefonprefix\usekomavar{durchwahl} \\}}%
	%\ifthenelse{\equal{\NTFDbrief@stil}{wagner}\OR\equal{\NTFDbrief@stil}{ernst}}{%
		%\ifkomavarempty{fromemail}{}{%
			%\l@addto@macro\NTFDlocation{\usekomavar*{fromemail} & \usekomavar{fromemail} \\}}%
	%}{}
	\l@addto@macro\NTFDrefnames{\\}%
	%\ifkomavarempty{datum}{%
		%\l@addto@macro\NTFDlocation{\usekomavar*{datum} & \today \\}}{%
		%\l@addto@macro\NTFDlocation{\usekomavar*{datum} & \usekomavar{datum} \\}}%
	\scriptsize\setlength{\unitlength}{1mm}%
	%\begin{picture}(0.1,0.1)(-50.3,0)
		%\put(0,0){\makebox(0,0)[br]{%
			%\renewcommand{\arraystretch}{1.5}%
			\begin{textblock*}{\useplength{refwidth}}(\useplength{refhpos},\useplength{refvpos})
		\begin{tabular}[b]{llll}%
			\NTFDrefnames
		\end{tabular}%
	\end{textblock*}
		%}}%
	%\end{picture}%
}}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brieffuß erste Seite
%
% Längen setzen
%\newlength{\NTFDbrief@fuss}
%\newlength{\NTFDbrief@kastenbreite}
%\ifthenelse{\equal{\NTFDbrief@stil}{intern}}{%
	%\setlength{\NTFDbrief@kastenbreite}{\TUBAFbausteineRandRi\TUBAFbausteine@bausteineeinheit}%
   %}{%
	%\setlength{\NTFDbrief@kastenbreite}{106.3988021pt}%
%}%
%\setlength{\NTFDbrief@fuss}{\paperwidth}%
%\addtolength{\NTFDbrief@fuss}{-\TUBAFbausteineRandLi\TUBAFbausteine@bausteineeinheit}% linker Rand
%\addtolength{\NTFDbrief@fuss}{-\NTFDbrief@kastenbreite}	% Breite Identkasten
%\ifthenelse{\equal{\NTFDbrief@stil}{intern}}{}{
	%\addtolength{\NTFDbrief@fuss}{-1.5ex}			% Korrektur
%}%
%%
%\newcommand{\NTFD@fusstrenner}{~|~}
\setkomavar{firstfoot}{%
		%% Breite der Fuszeile
		%\setlength{\NTFD@temp}{182.483mm}% Korrektur, sonst ovefull hbox
		%\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
		%\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge
		%\FPround\NTFDbausteine@breite\NTFDbausteine@tmpa5
	%	
		%% Y position
		%\setlength{\NTFD@temp}{\paperheight}% Korrektur, sonst ovefull hbox
		%\FPset{\NTFDbausteine@basislaenge}{\number\NTFD@temp}
		%\FPdiv\NTFDbausteine@tmpa\NTFDbausteine@basislaenge\TUBAFbausteine@basislaenge
		%\FPround\NTFDbausteine@tmpa\NTFDbausteine@tmpa5
		%\FPset{\myY}{15.4}
		%\FPsub\NTFDbausteine@tmpa\NTFDbausteine@tmpa\myY
		%\setplength{firstfootvpos}{}<++>	
		%\begin{textblock}{\NTFDbausteine@breite}(\NTFD@xoffset,\NTFDbausteine@tmpa)
		%\begin{textblock}{\useplength{firstfootwidth}}(\useplength{firstfoothpos},\useplength{firstfootvpos})
			\renewcommand{\NTFDFootFG}{NTFDgrau}
			\begin{tikzpicture}
				\NTFDFusszeile
			\end{tikzpicture}
		%\end{textblock}

	%\ifthenelse{\equal{\NTFDbrief@stil}{wagner}\OR\equal{\NTFDbrief@stil}{ernst}}{}{%
	%\hspace*{\TUBAFbausteineRandLi\TUBAFbausteine@bausteineeinheit}%
	%\begin{minipage}[b]{\NTFDbrief@fuss}
		%\fontsize{7}{9}\selectfont
		%\flushright%
		%~\hfill%
		%\ifkomavarempty{fromaddress}{}{{\def\\{\NTFD@fusstrenner}\usekomavar{fromaddress}}%
			%\ifkomavarempty{fromphone}{\ifkomavarempty{fromfax}{}{\NTFD@fusstrenner}}{\NTFD@fusstrenner}}%
		%\ifkomavarempty{fromphone}{}{\usekomavar*{fromphone} \NTFDTelefonprefix\usekomavar{fromphone}%
			%\ifkomavarempty{fromfax}{}{\NTFD@fusstrenner}}%
		%\ifkomavarempty{fromfax}{}{\usekomavar*{fromfax} \NTFDTelefonprefix\usekomavar{fromfax}}%
		%\vphantom{()}%
		%\newline
		%\ifkomavarempty{fromemail}{}{\usekomavar*{fromemail}\usekomavar{fromemail}%
			 %\ifkomavarempty{internet}{}{\NTFD@fusstrenner}}%
		%\ifkomavarempty{internet}{}{\usekomavar*{internet} \NTFDInternetprefix\usekomavar{internet}}%
	%\end{minipage}\hfill~
}%}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Briefkopf Folgeseiten
\pagestyle{scrheadings} %
\setkomavar{nexthead}{%
	\vspace*{0pt}\begin{minipage}[t]{\linewidth}%
		\footnotesize%
		\normalfont%
		\ifkomavarempty{subject}{}{\usekomavar{subject}\newline}
		\usekomavar{fromname} \hfill 
		\ifkomavarempty{datum}{\today}{\usekomavar{datum}} \hfill	
		\pagename\pageseparator\thepage\usekomavar*{interpage}\pagesatthisletter%
		\vskip0pt
		\ifthenelse{\isundefined{\NTFDhausschrift@strichstaerke}}{%
			\rule[1.5ex]{\linewidth}{.08ex}%
		   }{%
			\rule[1.5ex]{\linewidth}{.66\NTFDhausschrift@strichstaerke}%
		}
		\ifkomavarempty{subject}{}{\vspace*{-3ex}}% beseitigt overfull vbox
	\end{minipage}%
}%
%
\cfoot{%
	%\ifthenelse{\equal{\NTFDbrief@stil}{epost}}{%
		%~\hfill\raisebox{-15mm}{\NTFDkastenotb[fg=weiss]}\hspace*{-22.4mm}% elegantere Loesung? \NTFDkasten fnkt wg. textbox nicht
 	   %}{%
		%\ifthenelse{\isundefined{\NTFDhausschrift@strichstaerke}}{%
			%\rule{\linewidth}{.08ex}%
		   %}{%
			%\rule{\linewidth}{.66\NTFDhausschrift@strichstaerke}%
		%}%
	%}%
	% Breite der Fuszeile
	\rule[0pt]{0pt}{15.7pt}
	\renewcommand{\NTFDFootFG}{NTFDgrau}
	\begin{tikzpicture}
		\NTFDFusszeile
	\end{tikzpicture}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geschaeftszeile und Datum leeren
\removereffields%	löscht die horizontale Konstruktion über der Betreffzeile
\date{}%Datum leeren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backaddressseparator
\setkomavar{backaddressseparator}{~$\cdot$~}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die Unterschrift (signatur}
\newlength{\NTFDbrief@sigindent}
%
\setkomavar{signature}{%
	\ifkomavarempty{funktion}{%
		\usekomavar{fromname}%
	   }{%
		\ifthenelse{\equal{\NTFDbrief@funktion}{beides}\OR\equal{\NTFDbrief@funktion}{unten}}{%
			\settowidth{\@tempdima}{\usekomavar{fromname}}%
			\settowidth{\@tempdimb}{\small\usekomavar{funktion}}%
			\ifdim\@tempdima>\@tempdimb%
				\setlength{\@tempdimc}{\the\@tempdima}%
			   \else%
				\setlength{\@tempdimc}{\the\@tempdimb}%
				\setlength{\NTFDbrief@sigindent}{\the\@tempdimb}%
				\addtolength{\NTFDbrief@sigindent}{-\the\@tempdima}%
			\fi%
			\parbox{\@tempdimc}{%
				\centering%
				\hskip.25\NTFDbrief@sigindent\usekomavar{fromname}\newline%
				{\small\usekomavar{funktion}}%
			}%
		   }{%
			\usekomavar{fromname}
		}%
	}%
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neue Variablen
%
% unserschreiben
\newkomavar[%
	\usekomavar*{unserschreibenseparator}%
	\usekomavar{unserschreibenseparator}%
	]{unserschreiben}
\newkomavar[\unserschreibenname]{unserschreibenseparator}
\setkomavar{unserschreibenseparator}{:}
\newcommand{\unserschreibenname}{Unser Schreiben vom}
%
% datum
\newkomavar[%
	\usekomavar*{datumseparator}%
	\usekomavar{datumseparator}%
	]{datum}
\newkomavar[\datumname]{datumseparator}
\setkomavar{datumseparator}{:}
\newcommand{\datumname}{Datum}
%
% durchwahl
\newkomavar[%
	\usekomavar*{durchwahlseparator}%
	\usekomavar{durchwahlseparator}%
	]{durchwahl}
\newkomavar[\durchwahlname]{durchwahlseparator}
\setkomavar{durchwahlseparator}{:}
\newcommand{\durchwahlname}{Durchwahl}
%
% Internet
\newkomavar[%
	\usekomavar*{internetseparator}%
	\usekomavar{internetseparator}%
	]{internet}
\newkomavar[\internetname]{internetseparator}
\setkomavar{internetseparator}{:}
\newcommand{\internetname}{Internet}
%
% lehrstuhl
\newkomavar{lehrstuhl}
%
% institut
\newkomavar{institut}
%
% fakulaet
\newkomavar{fakultaet}
%
% funktion
\newkomavar{funktion}
%
% Seitenzahlzwischenwort
\newkomavar[%
	\usekomavar*{interpageseparator}%
	\usekomavar{interpageseparator}%
	]{interpage}
\newkomavar[\interpagename]{interpageseparator}
\setkomavar{interpageseparator}{~}
\newcommand{\interpagename}{\usekomavar{interpageseparator}von}
\newcommand{\pageseparator}{~}
%
% Postscriptum
\newkomavar[%
	\usekomavar*{psseparator}%
	\usekomavar{psseparator}%
	]{ps}
\newkomavar[\psname]{psseparator}
\setkomavar{psseparator}{:~}
\newcommand{\psname}{PS}
% Umdefinieren von \ps
\renewcommand*{\ps}[1]{\par%
	  \ifdim\parskip=\z@%
	    \vskip\baselineskip%
	  \fi\noindent%
	  \begingroup
	    \parbox[t]{\textwidth}{%
	      \ifkomavarempty{psseparator}{}{%
	        \@hangfrom{%
	          \strut\usekomavar*{psseparator}\usekomavar{psseparator}}%
	        }%
	      \ignorespaces #1\strut}%
	    \setlength{\parfillskip}{\z@ \@plus 1fil}\par
	  \endgroup
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ausrichtung Grussformel und Signatur: linksbündig
\renewcommand*{\raggedsignature}{\raggedright}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% für die Benutzerfreundlichkeit werden alle Variablen in die \NTFD... Form überführt
%
\newcommand{\NTFDBetreff}[1]{\setkomavar{subject}{#1}}
\newcommand{\NTFDIhrSchreiben}[1]{\setkomavar{yourmail}{#1}}
\newcommand{\NTFDIhrZeichen}[1]{\setkomavar{yourref}{#1}}
\NTFDIhrZeichen{\null}
\newcommand{\NTFDUnserSchreiben}[1]{\setkomavar{unserschreiben}{#1}}
\newcommand{\NTFDUnserZeichen}[1]{\setkomavar{myref}{#1}}
\newcommand{\NTFDDatum}[1]{\setkomavar{datum}{#1}}
\newcommand{\NTFDAbsendername}[1]{\setkomavar{fromname}{#1}}
	\AtEndLetter{\ifkomavarempty{fromname}{%
		\ClassError{\NTFDbrief@name}{Please fill \protect\NTFDAbsendername\protect{\protect}!}{}}{}}
\newcommand{\NTFDAbsenderanschrift}[1]{\setkomavar{fromaddress}{#1}}
\newcommand{\NTFDAbsendertelefon}[1]{\setkomavar{fromphone}{#1}}
\newcommand{\NTFDAbsenderfax}[1]{\setkomavar{fromfax}{#1}}
\newcommand{\NTFDAbsenderdurchwahl}[1]{\setkomavar{durchwahl}{#1}}
\newcommand{\NTFDEmail}[1]{\setkomavar{fromemail}{#1}}
\newcommand{\NTFDInternet}[1]{\setkomavar{internet}{#1}}
\newcommand{\NTFDFakultaet}[1]{\setkomavar{fakultaet}{#1}}
\newcommand{\NTFDInstitut}[1]{\setkomavar{institut}{#1}}
\newcommand{\NTFDLehrstuhl}[1]{\setkomavar{lehrstuhl}{#1}}
\newcommand{\NTFDFunktion}[1]{\setkomavar{funktion}{#1}}
%
\newcommand{\NTFDAbsenderdaten}[1]{\LoadLetterOption{#1}}
\newcommand{\NTFDAnrede}[1]{\opening{#1}}
\newcommand{\NTFDGrussformel}[1]{\ifthenelse{\equal{#1}{}}{\vspace*{3ex}\closing{Gl\"uckauf!}}{\vspace*{3ex}\closing{#1}}}
%
\newcommand{\NTFDVersandart}[1]{\setkomavar{specialmail}{#1}}
\newcommand{\NTFDPostscriptum}[1]{\ps{#1}}
\newcommand{\NTFDAnlagen}[1]{\encl{#1}}
\newcommand{\NTFDVerteiler}[1]{\cc{#1}}
%
\newcommand{\NTFDTelefonprefix}{+49--(0)3731--39~}
\newcommand{\NTFDInternetprefix}{http:/$\!$/}
%
\newcommand{\NTFDRetourAdresse}[1]{%
	\ifthenelse{\equal{#1}{}}{%
		\ifthenelse{\equal{\NTFDbrief@stil}{wagner}}{%
			\setkomavar{backaddress}{}%
		   }{%
			\setkomavar{backaddress}{TU Bergakademie Freiberg\usekomavar{backaddressseparator}%
				Akademiestr. 6\usekomavar{backaddressseparator}%
				D-09596 Freiberg}%
		}%
	   }{%
		\ifthenelse{\equal{\NTFDbrief@stil}{wagner}}{%
			\setkomavar{backaddress}{}%
		   }{%
	   		\setkomavar{backaddress}{#1}%
	   	}%
	}%
}
% wird nochmal separat gesetzt, falls kein \NTFDRetourAdresse angegeben
\ifthenelse{\equal{\NTFDbrief@stil}{wagner}}{\setkomavar{backaddress}{}}{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\endinput
