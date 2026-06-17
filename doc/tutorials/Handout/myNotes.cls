\NeedsTeXFormat{LaTeX2e} \ProvidesClass{myNotes}[10/06/2005 class for my notes]

\RequirePackage{setspace}
\newif\ifsinglesp@
\newif\ifdoublesp@
\doublesp@false

\DeclareOption{singlespace}{\singlesp@true}
\DeclareOption{doublespace}{\doublesp@true}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions

\LoadClass[a4paper]{article}

\hoffset0cm
\oddsidemargin0cm
\evensidemargin0cm
\textwidth16.7cm

% \oddsidemargin 8.6truemm \evensidemargin 0truemm
% \marginparwidth 40pt \marginparsep 10pt
% \topmargin 0truemm \headsep 10truemm
% \textheight 224truemm \textwidth 150truemm

\brokenpenalty=10000

\AtBeginDocument{%
\ifdoublesp@%
  \doublespacing%
\else\ifsinglesp@%
  \singlespacing%
\else%
  \onehalfspacing%
\fi%
\fi
}

\def\date#1{\gdef\@date{#1}}
\def\town#1{\gdef\@town{#1}}
\def\abstract#1{\gdef\@abstract{#1}}

\def\theequation{\arabic{section}.\arabic{equation}}

\newcommand {\Zb} {\mathbb{Z}} % integer numbers
\newcommand {\Qb} {\mathbb{Q}} % rartional numbers
\newcommand {\Rb} {\mathbb{R}} % real numbers
\newcommand {\Cb} {\mathbb{C}} % complex numbers
\newcommand {\Nb} {\mathbb{N}} % natural numbers
\newcommand {\Db} {\mathbb{D}}
\newcommand {\Eb} {\mathbb{E}}
\newcommand {\Fb} {\mathbb{F}}
\newcommand {\Gb} {\mathbb{G}}
\newcommand {\Hb} {\mathbb{H}}
\newcommand {\Jb} {\mathbb{J}}
\newcommand {\Kb} {\mathbb{K}}
\newcommand {\Lb} {\mathbb{L}}
\newcommand {\const} {\operatorname{const.}}
\newcommand {\tr} {\operatorname{tr}}
\newcommand {\ord} {\mathcal{O}}
\newcommand {\Lag} {\mathcal{L}}
\newcommand {\QCD} {\operatorname{QCD}}
\newcommand {\BRS} {\operatorname{BRS}}
\newcommand {\YM} {\operatorname{YM}}
\newcommand {\kin} {\operatorname{kin}}
\newcommand {\inv} {\operatorname{inv}}
\newcommand {\GF} {\operatorname{GF}}
\newcommand {\FP} {\operatorname{FP}}
\newcommand {\Res}[1] {\operatorname{Res} \Big|_{#1}}
\newcommand {\I} {\operatorname{I}}
\newcommand {\II} {\operatorname{II}}
\newcommand {\III} {\operatorname{III}}
\newcommand {\id} {\mathbf{1}} %identity
% collides with menukeys
%\newcommand {\shift} {\hspace{3mm}}
\newcommand {\nshift} {\hspace{-3mm}}
\newcommand {\rct} {\tilde{r}_c}
\newcommand {\mirror} {\operatorname{mirror}}
\newcommand {\vp} {\varphi}

\def\titlepage{
    \thispagestyle{empty}
    \begin{center}
		\includegraphics[width=0.7\textwidth]{assets/OpenMS_5cm_300ppi_transparent.png}\\[0.5cm]
		\LARGE\textbf{User Tutorial}\\[5cm]		
		\vfill
		%\vspace{2cm}
                \begin{figure}
		\begin{minipage}{.5\textwidth}
			\centering
			\includegraphics[width=.8\linewidth, height=2cm]{assets/deNBI_Logo_rgb.png}
		\end{minipage}%
		\begin{minipage}{.5\textwidth}
			\centering
			\figstrut{5cm}{2cm}
			\includegraphics[width=.8\linewidth, height=5cm]{assets/elixir_germany_white_background.png}
		\end{minipage}
                \end{figure}	
		%\hfill\vspace{-2cm}\includegraphics[width=0.3\textwidth]{./graphics/assets/deNBI_Logo_rgb.png}\hfill\hspace{5cm}\vspace{1cm}\includegraphics[width=0.3\textwidth]{./graphics/assets/elixir_germany_white_background.png}\\
	\end{center}
	\newpage
    \thispagestyle{empty}
	\vfill
	\begin{center}
		\Large The OpenMS Developers\\[1cm]
%		presented by: \\
%		Julianus Pfeuffer, Timo Sachsenberg \\[2cm]
		% either put here the presenting authors for workshops or all authors
%       \textsc{\Large License}\\[1cm]
		\vfill
		\today 
		~\\ %hard space
		\textbf{Creative Commons Attribution 4.0 International}
		~\\ %hard space
		\textbf{(CC BY 4.0)}\\
	\end{center}
	\newpage
} % end titlepage

\def\firstpages{
    \pagenumbering{arabic}
    \pagestyle{plain}
    \titlepage
    \tableofcontents
	\newpage
	\pagestyle{plain}
}
