function out1 = levelquartic(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,x1,x2)
%LEVELQUARTIC
%    OUT1 = LEVELQUARTIC(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,X1,X2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    31-Jan-2017 13:11:30

out1 = -x1.*(-x1.*(-x1.*(A1.*B2.^4+A5.*B1.^4.*4.0-B1.^4.*x1-A2.*B1.*B2.^3-A1.*B1.^3.*B5.*4.0-A2.*B1.^3.*B4.*3.0-A3.*B1.^3.*B3.*2.0-A4.*B1.^3.*B2+A1.*B1.^3.*x2.*4.0+A1.*B1.^2.*B3.^2.*2.0+A3.*B1.^2.*B2.^2-A1.*B1.*B2.^2.*B3.*4.0+A1.*B1.^2.*B2.*B4.*4.0+A2.*B1.^2.*B2.*B3.*3.0)+A1.^2.*B3.^4+A5.^2.*B1.^4.*6.0-A3.^2.*B1.^3.*x2.*2.0+A1.^2.*B1.^2.*B5.^2.*6.0+A1.^2.*B2.^2.*B4.^2.*2.0+A2.^2.*B1.^2.*B4.^2.*3.0+A3.^2.*B1.^2.*B3.^2+A1.^2.*B1.^2.*x2.^2.*6.0+A1.*A5.*B2.^4.*3.0+A2.^2.*B1.*B3.^3+A4.^2.*B1.^3.*B3+A3.^2.*B1.^3.*B5.*2.0-A1.*A3.*B1.^2.*B4.^2.*3.0+A1.*A3.*B2.^2.*B3.^2+A1.*A5.*B1.^2.*B3.^2.*6.0-A2.*A4.*B1.^2.*B3.^2.*2.0+A3.*A5.*B1.^2.*B2.^2.*3.0+A1.^2.*B1.*B3.*B4.^2.*4.0-A1.^2.*B1.*B3.^2.*B5.*4.0-A1.^2.*B2.*B3.^2.*B4.*4.0+A2.^2.*B1.*B2.^2.*B5.*3.0-A3.^2.*B1.^2.*B2.*B4.*2.0+A1.^2.*B2.^2.*B3.*B5.*4.0-A2.^2.*B1.^2.*B3.*B5.*3.0+A1.^2.*B1.*B3.^2.*x2.*4.0-A2.^2.*B1.*B2.^2.*x2.*3.0-A1.^2.*B2.^2.*B3.*x2.*4.0+A2.^2.*B1.^2.*B3.*x2.*3.0-A1.^2.*B1.^2.*B5.*x2.*1.2e1-A1.*A2.*B2.*B3.^3-A1.*A3.*B1.*B3.^3.*2.0-A1.*A2.*B2.^3.*B5.*3.0-A1.*A3.*B2.^3.*B4.*2.0-A1.*A4.*B2.^3.*B3-A2.*A5.*B1.*B2.^3.*3.0-A1.*A5.*B1.^3.*B5.*1.2e1+A2.*A4.*B1.^3.*B5.*4.0-A2.*A5.*B1.^3.*B4.*9.0+A3.*A4.*B1.^3.*B4.*3.0-A3.*A5.*B1.^3.*B3.*6.0-A4.*A5.*B1.^3.*B2.*3.0+A1.*A2.*B2.^3.*x2.*3.0+A1.*A5.*B1.^3.*x2.*1.2e1-A2.*A4.*B1.^3.*x2.*4.0-A1.*A2.*B1.*B2.*B4.^2.*5.0+A1.*A2.*B1.*B3.^2.*B4+A1.*A4.*B1.*B2.*B3.^2.*3.0-A2.*A3.*B1.*B2.*B3.^2+A1.*A2.*B2.^2.*B3.*B4.*3.0+A1.*A3.*B1.*B2.^2.*B5.*2.0+A1.*A4.*B1.*B2.^2.*B4-A1.*A5.*B1.*B2.^2.*B3.*1.2e1+A2.*A3.*B1.*B2.^2.*B4.*2.0+A2.*A4.*B1.*B2.^2.*B3+A1.*A2.*B1.^2.*B4.*B5.*5.0+A1.*A3.*B1.^2.*B3.*B5.*2.0-A1.*A4.*B1.^2.*B2.*B5-A1.*A4.*B1.^2.*B3.*B4.*5.0+A1.*A5.*B1.^2.*B2.*B4.*1.2e1-A2.*A3.*B1.^2.*B2.*B5.*5.0+A2.*A3.*B1.^2.*B3.*B4-A2.*A4.*B1.^2.*B2.*B4+A2.*A5.*B1.^2.*B2.*B3.*9.0-A3.*A4.*B1.^2.*B2.*B3-A2.^2.*B1.*B2.*B3.*B4.*3.0-A1.^2.*B1.*B2.*B4.*B5.*8.0-A1.*A3.*B1.*B2.^2.*x2.*2.0-A1.*A2.*B1.^2.*B4.*x2.*5.0-A1.*A3.*B1.^2.*B3.*x2.*2.0+A1.*A4.*B1.^2.*B2.*x2+A2.*A3.*B1.^2.*B2.*x2.*5.0+A1.^2.*B1.*B2.*B4.*x2.*8.0-A1.*A2.*B1.*B2.*B3.*x2.*2.0+A1.*A2.*B1.*B2.*B3.*B5.*2.0+A1.*A3.*B1.*B2.*B3.*B4.*4.0)+A1.^3.*B4.^4+A5.^3.*B1.^4.*4.0+A1.^3.*B1.*x2.^3.*4.0+A3.^3.*B1.^2.*B4.^2+A1.^3.*B3.^2.*B5.^2.*2.0+A1.^3.*B3.^2.*x2.^2.*2.0+A1.*A5.^2.*B2.^4.*3.0+A1.^2.*A5.*B3.^4.*2.0-A1.^3.*B1.*B5.^3.*4.0-A2.^3.*B1.*B4.^3-A4.^3.*B1.^3.*B4+A1.*A2.^2.*B2.*B4.^3-A2.*A5.^2.*B1.*B2.^3.*3.0-A1.^2.*A2.*B3.*B4.^3-A1.^2.*A3.*B2.*B4.^3.*2.0-A1.^2.*A4.*B1.*B4.^3.*3.0+A1.*A4.^2.*B2.^3.*B4+A2.^2.*A5.*B1.*B3.^3.*2.0-A1.*A5.^2.*B1.^3.*B5.*1.2e1-A2.*A5.^2.*B1.^3.*B4.*9.0-A3.*A5.^2.*B1.^3.*B3.*6.0-A4.*A5.^2.*B1.^3.*B2.*3.0-A1.^2.*A3.*B3.^3.*B5.*2.0-A1.^2.*A4.*B3.^3.*B4-A3.*A4.^2.*B1.^3.*B5.*4.0+A4.^2.*A5.*B1.^3.*B3.*2.0+A3.^2.*A5.*B1.^3.*B5.*4.0-A2.^3.*B1.*B2.*B5.^2.*3.0+A1.^3.*B2.*B4.*B5.^2.*4.0-A3.^3.*B1.^2.*B3.*B5.*2.0-A1.^3.*B3.*B4.^2.*B5.*4.0+A1.*A5.^2.*B1.^3.*x2.*1.2e1+A1.^2.*A3.*B3.^3.*x2.*2.0+A3.*A4.^2.*B1.^3.*x2.*4.0-A3.^2.*A5.*B1.^3.*x2.*4.0-A2.^3.*B1.*B2.*x2.^2.*3.0-A1.^3.*B1.*B5.*x2.^2.*1.2e1+A1.^3.*B1.*B5.^2.*x2.*1.2e1+A1.^3.*B2.*B4.*x2.^2.*4.0+A3.^3.*B1.^2.*B3.*x2.*2.0+A1.^3.*B3.*B4.^2.*x2.*4.0-A1.^3.*B3.^2.*B5.*x2.*4.0+A1.*A2.^2.*B2.^2.*B5.^2.*3.0-A1.*A3.^2.*B1.^2.*B5.^2.*4.0+A1.*A3.^2.*B2.^2.*B4.^2+A1.*A4.^2.*B1.^2.*B4.^2.*3.0+A1.*A5.^2.*B1.^2.*B3.^2.*6.0+A3.*A5.^2.*B1.^2.*B2.^2.*3.0-A1.^2.*A3.*B2.^2.*B5.^2.*3.0+A1.^2.*A3.*B3.^2.*B4.^2+A2.^2.*A3.*B1.^2.*B5.^2.*4.0+A1.^2.*A5.*B1.^2.*B5.^2.*1.2e1+A1.^2.*A5.*B2.^2.*B4.^2.*4.0+A2.^2.*A5.*B1.^2.*B4.^2.*6.0+A3.^2.*A5.*B1.^2.*B3.^2.*2.0+A1.*A2.^2.*B2.^2.*x2.^2.*3.0-A1.*A3.^2.*B1.^2.*x2.^2.*4.0-A1.^2.*A3.*B2.^2.*x2.^2.*3.0+A2.^2.*A3.*B1.^2.*x2.^2.*4.0+A1.^2.*A5.*B1.^2.*x2.^2.*1.2e1+A1.*A2.*A3.*B1.*B4.^3.*3.0-A1.*A2.*A5.*B2.*B3.^3.*2.0-A1.*A3.*A5.*B1.*B3.^3.*4.0-A1.*A2.*A5.*B2.^3.*B5.*6.0+A1.*A3.*A4.*B2.^3.*B5.*3.0-A1.*A3.*A5.*B2.^3.*B4.*4.0-A1.*A4.*A5.*B2.^3.*B3.*2.0+A2.*A4.*A5.*B1.^3.*B5.*8.0+A3.*A4.*A5.*B1.^3.*B4.*6.0+A2.^3.*B1.*B3.*B4.*B5.*3.0+A1.*A2.*A5.*B2.^3.*x2.*6.0-A1.*A3.*A4.*B2.^3.*x2.*3.0-A2.*A4.*A5.*B1.^3.*x2.*8.0+A2.^3.*B1.*B2.*B5.*x2.*6.0-A2.^3.*B1.*B3.*B4.*x2.*3.0-A1.^3.*B2.*B4.*B5.*x2.*8.0-A1.*A2.*A4.*B1.^2.*B5.^2.*8.0-A1.*A2.*A4.*B2.^2.*B4.^2.*2.0-A1.*A3.*A5.*B1.^2.*B4.^2.*6.0+A1.*A3.*A5.*B2.^2.*B3.^2.*2.0-A2.*A3.*A4.*B1.^2.*B4.^2.*3.0-A2.*A4.*A5.*B1.^2.*B3.^2.*4.0+A1.*A2.^2.*B1.*B3.*B5.^2.*2.0-A1.*A3.^2.*B1.*B3.*B4.^2.*2.0-A1.*A5.^2.*B1.*B2.^2.*B3.*1.2e1-A2.*A3.^2.*B1.*B2.*B4.^2-A1.*A2.^2.*B1.*B4.^2.*B5+A1.*A3.^2.*B1.*B3.^2.*B5.*4.0-A1.*A4.^2.*B1.*B2.^2.*B5+A1.*A5.^2.*B1.^2.*B2.*B4.*1.2e1-A2.*A4.^2.*B1.*B2.^2.*B4+A2.*A5.^2.*B1.^2.*B2.*B3.*9.0-A1.^2.*A2.*B1.*B4.*B5.^2-A1.^2.*A2.*B2.*B3.*B5.^2.*5.0+A1.^2.*A3.*B1.*B3.*B5.^2.*2.0+A1.^2.*A4.*B1.*B2.*B5.^2.*5.0+A2.^2.*A3.*B1.*B3.*B4.^2+A2.^2.*A4.*B1.*B2.*B4.^2.*2.0-A1.*A3.^2.*B2.^2.*B3.*B5.*2.0+A1.*A4.^2.*B1.^2.*B3.*B5.*2.0+A2.*A4.^2.*B1.^2.*B2.*B5+A2.*A4.^2.*B1.^2.*B3.*B4.*2.0+A3.*A4.^2.*B1.^2.*B2.*B4+A1.^2.*A2.*B2.*B4.^2.*B5+A1.^2.*A3.*B1.*B4.^2.*B5.*2.0+A1.^2.*A4.*B2.*B3.*B4.^2.*3.0+A1.^2.*A5.*B1.*B3.*B4.^2.*8.0-A2.^2.*A3.*B1.*B3.^2.*B5.*2.0-A2.^2.*A4.*B1.*B3.^2.*B4+A2.*A3.^2.*B1.^2.*B4.*B5+A1.^2.*A2.*B3.^2.*B4.*B5.*3.0+A1.^2.*A4.*B2.*B3.^2.*B5-A1.^2.*A5.*B1.*B3.^2.*B5.*8.0-A1.^2.*A5.*B2.*B3.^2.*B4.*8.0+A2.^2.*A5.*B1.*B2.^2.*B5.*6.0+A3.^2.*A4.*B1.^2.*B2.*B5.*3.0-A3.^2.*A4.*B1.^2.*B3.*B4-A3.^2.*A5.*B1.^2.*B2.*B4.*4.0-A1.^2.*A4.*B2.^2.*B4.*B5.*5.0+A1.^2.*A5.*B2.^2.*B3.*B5.*8.0-A2.^2.*A4.*B1.^2.*B4.*B5.*5.0-A2.^2.*A5.*B1.^2.*B3.*B5.*6.0-A1.*A2.*A4.*B1.^2.*x2.^2.*8.0+A1.*A2.^2.*B1.*B3.*x2.^2.*2.0+A1.*A2.^2.*B1.*B4.^2.*x2-A1.*A3.^2.*B1.*B3.^2.*x2.*4.0+A1.*A4.^2.*B1.*B2.^2.*x2-A1.^2.*A2.*B1.*B4.*x2.^2-A1.^2.*A2.*B2.*B3.*x2.^2.*5.0+A1.^2.*A3.*B1.*B3.*x2.^2.*2.0+A1.^2.*A4.*B1.*B2.*x2.^2.*5.0+A1.*A3.^2.*B2.^2.*B3.*x2.*2.0-A1.*A4.^2.*B1.^2.*B3.*x2.*2.0-A2.*A4.^2.*B1.^2.*B2.*x2-A1.^2.*A2.*B2.*B4.^2.*x2-A1.^2.*A3.*B1.*B4.^2.*x2.*2.0+A2.^2.*A3.*B1.*B3.^2.*x2.*2.0-A1.*A2.^2.*B2.^2.*B5.*x2.*6.0+A1.*A3.^2.*B1.^2.*B5.*x2.*8.0-A2.*A3.^2.*B1.^2.*B4.*x2-A1.^2.*A2.*B3.^2.*B4.*x2.*3.0-A1.^2.*A4.*B2.*B3.^2.*x2+A1.^2.*A5.*B1.*B3.^2.*x2.*8.0-A2.^2.*A5.*B1.*B2.^2.*x2.*6.0-A3.^2.*A4.*B1.^2.*B2.*x2.*3.0+A1.^2.*A3.*B2.^2.*B5.*x2.*6.0+A1.^2.*A4.*B2.^2.*B4.*x2.*5.0-A1.^2.*A5.*B2.^2.*B3.*x2.*8.0-A2.^2.*A3.*B1.^2.*B5.*x2.*8.0+A2.^2.*A4.*B1.^2.*B4.*x2.*5.0+A2.^2.*A5.*B1.^2.*B3.*x2.*6.0-A1.^2.*A5.*B1.^2.*B5.*x2.*2.4e1+A1.*A2.*A3.*B1.*B2.*B5.^2.*2.0-A1.*A2.*A3.*B2.*B3.*B4.^2-A1.*A2.*A4.*B1.*B3.*B4.^2-A1.*A2.*A5.*B1.*B2.*B4.^2.*1.0e1+A1.*A3.*A4.*B1.*B2.*B4.^2+A1.*A2.*A3.*B2.*B3.^2.*B5.*2.0+A1.*A2.*A4.*B2.*B3.^2.*B4+A1.*A2.*A5.*B1.*B3.^2.*B4.*2.0+A1.*A3.*A4.*B1.*B3.^2.*B4.*2.0+A1.*A4.*A5.*B1.*B2.*B3.^2.*6.0-A2.*A3.*A5.*B1.*B2.*B3.^2.*2.0+A1.*A2.*A3.*B2.^2.*B4.*B5-A1.*A2.*A4.*B2.^2.*B3.*B5+A1.*A2.*A5.*B2.^2.*B3.*B4.*6.0-A1.*A3.*A4.*B2.^2.*B3.*B4+A1.*A3.*A5.*B1.*B2.^2.*B5.*4.0+A1.*A4.*A5.*B1.*B2.^2.*B4.*2.0-A2.*A3.*A4.*B1.*B2.^2.*B5.*3.0+A2.*A3.*A5.*B1.*B2.^2.*B4.*4.0+A2.*A4.*A5.*B1.*B2.^2.*B3.*2.0+A1.*A2.*A5.*B1.^2.*B4.*B5.*1.0e1+A1.*A3.*A4.*B1.^2.*B4.*B5.*2.0+A1.*A3.*A5.*B1.^2.*B3.*B5.*4.0-A1.*A4.*A5.*B1.^2.*B2.*B5.*2.0-A1.*A4.*A5.*B1.^2.*B3.*B4.*1.0e1+A2.*A3.*A4.*B1.^2.*B3.*B5.*4.0-A2.*A3.*A5.*B1.^2.*B2.*B5.*1.0e1+A2.*A3.*A5.*B1.^2.*B3.*B4.*2.0-A2.*A4.*A5.*B1.^2.*B2.*B4.*2.0-A3.*A4.*A5.*B1.^2.*B2.*B3.*2.0-A1.*A4.^2.*B1.*B2.*B3.*B4.*3.0+A2.*A3.^2.*B1.*B2.*B3.*B5.*2.0-A1.*A2.^2.*B2.*B3.*B4.*B5.*3.0-A2.^2.*A3.*B1.*B2.*B4.*B5+A2.^2.*A4.*B1.*B2.*B3.*B5-A2.^2.*A5.*B1.*B2.*B3.*B4.*6.0+A1.^2.*A3.*B2.*B3.*B4.*B5.*4.0+A1.^2.*A4.*B1.*B3.*B4.*B5.*2.0-A1.^2.*A5.*B1.*B2.*B4.*B5.*1.6e1+A1.*A2.*A3.*B1.*B2.*x2.^2.*2.0-A1.*A2.*A3.*B2.*B3.^2.*x2.*2.0-A1.*A2.*A3.*B2.^2.*B4.*x2+A1.*A2.*A4.*B2.^2.*B3.*x2-A1.*A3.*A5.*B1.*B2.^2.*x2.*4.0+A2.*A3.*A4.*B1.*B2.^2.*x2.*3.0+A1.*A2.*A4.*B1.^2.*B5.*x2.*1.6e1-A1.*A2.*A5.*B1.^2.*B4.*x2.*1.0e1-A1.*A3.*A4.*B1.^2.*B4.*x2.*2.0-A1.*A3.*A5.*B1.^2.*B3.*x2.*4.0+A1.*A4.*A5.*B1.^2.*B2.*x2.*2.0-A2.*A3.*A4.*B1.^2.*B3.*x2.*4.0+A2.*A3.*A5.*B1.^2.*B2.*x2.*1.0e1-A2.*A3.^2.*B1.*B2.*B3.*x2.*2.0-A1.*A2.^2.*B1.*B3.*B5.*x2.*4.0+A1.*A2.^2.*B2.*B3.*B4.*x2.*3.0+A2.^2.*A3.*B1.*B2.*B4.*x2-A2.^2.*A4.*B1.*B2.*B3.*x2+A1.^2.*A2.*B1.*B4.*B5.*x2.*2.0+A1.^2.*A2.*B2.*B3.*B5.*x2.*1.0e1-A1.^2.*A3.*B1.*B3.*B5.*x2.*4.0-A1.^2.*A3.*B2.*B3.*B4.*x2.*4.0-A1.^2.*A4.*B1.*B2.*B5.*x2.*1.0e1-A1.^2.*A4.*B1.*B3.*B4.*x2.*2.0+A1.^2.*A5.*B1.*B2.*B4.*x2.*1.6e1-A1.*A2.*A3.*B1.*B3.*B4.*B5.*8.0+A1.*A2.*A4.*B1.*B2.*B4.*B5.*1.0e1+A1.*A2.*A5.*B1.*B2.*B3.*B5.*4.0-A1.*A3.*A4.*B1.*B2.*B3.*B5.*8.0+A1.*A3.*A5.*B1.*B2.*B3.*B4.*8.0+A2.*A3.*A4.*B1.*B2.*B3.*B4-A1.*A2.*A3.*B1.*B2.*B5.*x2.*4.0+A1.*A2.*A3.*B1.*B3.*B4.*x2.*8.0-A1.*A2.*A4.*B1.*B2.*B4.*x2.*1.0e1-A1.*A2.*A5.*B1.*B2.*B3.*x2.*4.0+A1.*A3.*A4.*B1.*B2.*B3.*x2.*8.0)+A1.^4.*B5.^4+A5.^4.*B1.^4+A1.^4.*x2.^4-A2.^4.*B1.*x2.^3-A4.^4.*B1.^3.*x2-A1.^4.*B5.*x2.^3.*4.0-A1.^4.*B5.^3.*x2.*4.0+A1.^2.*A5.^2.*B3.^4+A3.^4.*B1.^2.*B5.^2+A3.^4.*B1.^2.*x2.^2+A1.^4.*B5.^2.*x2.^2.*6.0+A1.*A5.^3.*B2.^4+A1.^3.*A5.*B4.^4+A2.^4.*B1.*B5.^3+A4.^4.*B1.^3.*B5+A1.^2.*A3.^2.*B3.^2.*x2.^2+A1.^2.*A4.^2.*B2.^2.*x2.^2.*3.0+A1.^2.*A5.^2.*B1.^2.*x2.^2.*6.0+A2.^2.*A4.^2.*B1.^2.*x2.^2.*2.0-A1.*A2.^3.*B2.*B5.^3-A2.*A5.^3.*B1.*B2.^3-A1.*A4.^3.*B2.^3.*B5-A1.*A5.^3.*B1.^3.*B5.*4.0-A2.*A5.^3.*B1.^3.*B4.*3.0-A3.*A5.^3.*B1.^3.*B3.*2.0-A4.*A5.^3.*B1.^3.*B2-A1.^3.*A2.*B4.*B5.^3-A1.^3.*A3.*B3.*B5.^3.*2.0-A1.^3.*A4.*B2.*B5.^3.*3.0-A1.^3.*A5.*B1.*B5.^3.*4.0-A2.^3.*A5.*B1.*B4.^3-A1.^3.*A4.*B4.^3.*B5-A4.^3.*A5.*B1.^3.*B4+A1.*A2.^3.*B2.*x2.^3+A1.*A4.^3.*B2.^3.*x2+A1.*A5.^3.*B1.^3.*x2.*4.0+A1.^3.*A2.*B4.*x2.^3+A1.^3.*A3.*B3.*x2.^3.*2.0+A1.^3.*A4.*B2.*x2.^3.*3.0+A1.^3.*A5.*B1.*x2.^3.*4.0+A1.^3.*A4.*B4.^3.*x2+A2.^4.*B1.*B5.*x2.^2.*3.0-A2.^4.*B1.*B5.^2.*x2.*3.0-A3.^4.*B1.^2.*B5.*x2.*2.0+A1.*A5.^3.*B1.^2.*B3.^2.*2.0+A1.^2.*A3.^2.*B1.*B5.^3.*2.0+A1.*A3.^3.*B2.^2.*B5.^2+A3.*A5.^3.*B1.^2.*B2.^2+A1.^2.*A2.^2.*B3.*B5.^3+A2.^2.*A5.^2.*B1.*B3.^3+A1.^2.*A4.^2.*B3.^3.*B5+A1.^3.*A3.*B4.^2.*B5.^2+A3.^3.*A5.*B1.^2.*B4.^2+A4.^2.*A5.^2.*B1.^3.*B3+A1.^3.*A5.*B3.^2.*B5.^2.*2.0+A3.^2.*A5.^2.*B1.^3.*B5.*2.0-A1.^2.*A3.^2.*B1.*x2.^3.*2.0+A1.*A3.^3.*B2.^2.*x2.^2-A1.^2.*A2.^2.*B3.*x2.^3-A1.^2.*A4.^2.*B3.^3.*x2+A1.^3.*A3.*B4.^2.*x2.^2+A1.^3.*A5.*B3.^2.*x2.^2.*2.0-A3.^2.*A5.^2.*B1.^3.*x2.*2.0+A1.^2.*A3.^2.*B3.^2.*B5.^2+A1.^2.*A4.^2.*B2.^2.*B5.^2.*3.0+A1.^2.*A5.^2.*B1.^2.*B5.^2.*6.0+A1.^2.*A5.^2.*B2.^2.*B4.^2.*2.0+A2.^2.*A4.^2.*B1.^2.*B5.^2.*2.0+A2.^2.*A5.^2.*B1.^2.*B4.^2.*3.0+A3.^2.*A5.^2.*B1.^2.*B3.^2-A1.*A2.^2.*A3.*B1.*B5.^3.*4.0-A1.*A2.*A5.^2.*B2.*B3.^3-A1.*A3.*A5.^2.*B1.*B3.^3.*2.0+A1.^2.*A2.*A3.*B2.*B5.^3.*3.0+A1.^2.*A2.*A4.*B1.*B5.^3.*4.0+A1.*A2.^2.*A5.*B2.*B4.^3-A1.*A2.*A5.^2.*B2.^3.*B5.*3.0-A1.*A3.*A5.^2.*B2.^3.*B4.*2.0-A1.*A4.*A5.^2.*B2.^3.*B3-A1.^2.*A2.*A5.*B3.*B4.^3-A1.^2.*A3.*A5.*B2.*B4.^3.*2.0-A1.^2.*A4.*A5.*B1.*B4.^3.*3.0+A1.*A4.^2.*A5.*B2.^3.*B4+A2.*A4.*A5.^2.*B1.^3.*B5.*4.0+A3.*A4.*A5.^2.*B1.^3.*B4.*3.0-A1.^2.*A3.*A5.*B3.^3.*B5.*2.0-A1.^2.*A4.*A5.*B3.^3.*B4-A3.*A4.^2.*A5.*B1.^3.*B5.*4.0-A1.*A5.^3.*B1.*B2.^2.*B3.*4.0-A1.*A3.^3.*B1.*B3.*B5.^2.*2.0+A1.*A5.^3.*B1.^2.*B2.*B4.*4.0-A2.*A3.^3.*B1.*B2.*B5.^2+A2.*A5.^3.*B1.^2.*B2.*B3.*3.0+A2.*A4.^3.*B1.*B2.^2.*B5-A1.*A4.^3.*B1.^2.*B4.*B5.*3.0-A2.*A4.^3.*B1.^2.*B3.*B5.*2.0-A3.*A4.^3.*B1.^2.*B2.*B5-A2.^3.*A3.*B1.*B4.*B5.^2-A2.^3.*A4.*B1.*B3.*B5.^2.*2.0-A2.^3.*A5.*B1.*B2.*B5.^2.*3.0+A2.^3.*A4.*B1.*B4.^2.*B5+A1.^3.*A4.*B3.*B4.*B5.^2.*3.0+A1.^3.*A5.*B2.*B4.*B5.^2.*4.0-A3.^3.*A4.*B1.^2.*B4.*B5-A3.^3.*A5.*B1.^2.*B3.*B5.*2.0-A1.^3.*A5.*B3.*B4.^2.*B5.*4.0+A1.*A2.^2.*A3.*B1.*x2.^3.*4.0-A1.^2.*A2.*A3.*B2.*x2.^3.*3.0-A1.^2.*A2.*A4.*B1.*x2.^3.*4.0+A1.*A2.*A5.^2.*B2.^3.*x2.*3.0-A2.*A4.*A5.^2.*B1.^3.*x2.*4.0+A1.^2.*A3.*A5.*B3.^3.*x2.*2.0+A3.*A4.^2.*A5.*B1.^3.*x2.*4.0-A1.*A3.^3.*B1.*B3.*x2.^2.*2.0-A2.*A3.^3.*B1.*B2.*x2.^2-A2.*A4.^3.*B1.*B2.^2.*x2-A1.*A2.^3.*B2.*B5.*x2.^2.*3.0+A1.*A2.^3.*B2.*B5.^2.*x2.*3.0+A1.*A4.^3.*B1.^2.*B4.*x2.*3.0+A2.*A4.^3.*B1.^2.*B3.*x2.*2.0+A3.*A4.^3.*B1.^2.*B2.*x2-A2.^3.*A3.*B1.*B4.*x2.^2-A2.^3.*A4.*B1.*B3.*x2.^2.*2.0-A2.^3.*A5.*B1.*B2.*x2.^2.*3.0-A1.*A3.^3.*B2.^2.*B5.*x2.*2.0-A2.^3.*A4.*B1.*B4.^2.*x2-A1.^3.*A2.*B4.*B5.*x2.^2.*3.0+A1.^3.*A2.*B4.*B5.^2.*x2.*3.0-A1.^3.*A3.*B3.*B5.*x2.^2.*6.0+A1.^3.*A3.*B3.*B5.^2.*x2.*6.0-A1.^3.*A4.*B2.*B5.*x2.^2.*9.0+A1.^3.*A4.*B2.*B5.^2.*x2.*9.0+A1.^3.*A4.*B3.*B4.*x2.^2.*3.0-A1.^3.*A5.*B1.*B5.*x2.^2.*1.2e1+A1.^3.*A5.*B1.*B5.^2.*x2.*1.2e1+A1.^3.*A5.*B2.*B4.*x2.^2.*4.0+A3.^3.*A4.*B1.^2.*B4.*x2+A3.^3.*A5.*B1.^2.*B3.*x2.*2.0-A1.^3.*A3.*B4.^2.*B5.*x2.*2.0+A1.^3.*A5.*B3.*B4.^2.*x2.*4.0-A1.^3.*A5.*B3.^2.*B5.*x2.*4.0+A1.*A3.*A4.^2.*B1.^2.*B5.^2.*4.0-A1.*A3.*A5.^2.*B1.^2.*B4.^2.*3.0+A1.*A3.*A5.^2.*B2.^2.*B3.^2+A1.*A2.^2.*A5.*B2.^2.*B5.^2.*3.0-A1.*A3.^2.*A5.*B1.^2.*B5.^2.*4.0+A1.*A3.^2.*A5.*B2.^2.*B4.^2+A1.*A4.^2.*A5.*B1.^2.*B4.^2.*3.0-A2.*A4.*A5.^2.*B1.^2.*B3.^2.*2.0-A2.*A3.^2.*A4.*B1.^2.*B5.^2.*4.0-A1.^2.*A2.*A4.*B3.^2.*B5.^2.*2.0-A1.^2.*A3.*A5.*B2.^2.*B5.^2.*3.0+A1.^2.*A3.*A5.*B3.^2.*B4.^2+A2.^2.*A3.*A5.*B1.^2.*B5.^2.*4.0-A1.^2.*A4.^2.*B1.*B3.*B5.^2.*3.0+A1.^2.*A5.^2.*B1.*B3.*B4.^2.*4.0+A2.^2.*A3.^2.*B1.*B3.*B5.^2-A1.^2.*A3.^2.*B2.*B4.*B5.^2.*2.0+A1.^2.*A4.^2.*B1.*B4.^2.*B5.*3.0-A1.^2.*A5.^2.*B1.*B3.^2.*B5.*4.0-A1.^2.*A5.^2.*B2.*B3.^2.*B4.*4.0+A2.^2.*A4.^2.*B1.*B3.^2.*B5+A2.^2.*A5.^2.*B1.*B2.^2.*B5.*3.0-A3.^2.*A5.^2.*B1.^2.*B2.*B4.*2.0+A1.^2.*A5.^2.*B2.^2.*B3.*B5.*4.0-A2.^2.*A5.^2.*B1.^2.*B3.*B5.*3.0+A3.^2.*A4.^2.*B1.^2.*B3.*B5+A1.*A3.*A4.^2.*B1.^2.*x2.^2.*4.0+A1.*A2.^2.*A5.*B2.^2.*x2.^2.*3.0-A1.*A3.^2.*A5.*B1.^2.*x2.^2.*4.0-A2.*A3.^2.*A4.*B1.^2.*x2.^2.*4.0-A1.^2.*A2.*A4.*B3.^2.*x2.^2.*2.0-A1.^2.*A3.*A5.*B2.^2.*x2.^2.*3.0+A2.^2.*A3.*A5.*B1.^2.*x2.^2.*4.0-A1.^2.*A4.^2.*B1.*B3.*x2.^2.*3.0+A2.^2.*A3.^2.*B1.*B3.*x2.^2+A1.^2.*A3.^2.*B1.*B5.*x2.^2.*6.0-A1.^2.*A3.^2.*B1.*B5.^2.*x2.*6.0-A1.^2.*A3.^2.*B2.*B4.*x2.^2.*2.0-A1.^2.*A4.^2.*B1.*B4.^2.*x2.*3.0+A1.^2.*A5.^2.*B1.*B3.^2.*x2.*4.0-A2.^2.*A4.^2.*B1.*B3.^2.*x2-A2.^2.*A5.^2.*B1.*B2.^2.*x2.*3.0+A1.^2.*A2.^2.*B3.*B5.*x2.^2.*3.0-A1.^2.*A2.^2.*B3.*B5.^2.*x2.*3.0-A1.^2.*A5.^2.*B2.^2.*B3.*x2.*4.0+A2.^2.*A5.^2.*B1.^2.*B3.*x2.*3.0-A3.^2.*A4.^2.*B1.^2.*B3.*x2-A1.^2.*A3.^2.*B3.^2.*B5.*x2.*2.0-A1.^2.*A4.^2.*B2.^2.*B5.*x2.*6.0-A1.^2.*A5.^2.*B1.^2.*B5.*x2.*1.2e1-A2.^2.*A4.^2.*B1.^2.*B5.*x2.*4.0+A1.*A2.*A3.*A5.*B1.*B4.^3.*3.0+A1.*A3.*A4.*A5.*B2.^3.*B5.*3.0+A1.*A4.^3.*B1.*B2.*B3.*B5.*3.0+A2.^3.*A5.*B1.*B3.*B4.*B5.*3.0-A1.*A3.*A4.*A5.*B2.^3.*x2.*3.0-A1.*A4.^3.*B1.*B2.*B3.*x2.*3.0+A1.*A3.^3.*B1.*B3.*B5.*x2.*4.0+A2.*A3.^3.*B1.*B2.*B5.*x2.*2.0+A2.^3.*A3.*B1.*B4.*B5.*x2.*2.0+A2.^3.*A4.*B1.*B3.*B5.*x2.*4.0+A2.^3.*A5.*B1.*B2.*B5.*x2.*6.0-A2.^3.*A5.*B1.*B3.*B4.*x2.*3.0-A1.^3.*A4.*B3.*B4.*B5.*x2.*6.0-A1.^3.*A5.*B2.*B4.*B5.*x2.*8.0-A1.*A2.*A3.*A4.*B2.^2.*B5.^2.*3.0-A1.*A2.*A4.*A5.*B1.^2.*B5.^2.*8.0-A1.*A2.*A4.*A5.*B2.^2.*B4.^2.*2.0-A2.*A3.*A4.*A5.*B1.^2.*B4.^2.*3.0-A1.*A2.*A4.^2.*B1.*B2.*B5.^2.*5.0-A1.*A2.*A5.^2.*B1.*B2.*B4.^2.*5.0+A1.*A2.*A3.^2.*B1.*B4.*B5.^2.*3.0-A1.*A2.*A3.^2.*B2.*B3.*B5.^2+A1.*A2.*A5.^2.*B1.*B3.^2.*B4+A1.*A4.*A5.^2.*B1.*B2.*B3.^2.*3.0+A1.*A3.^2.*A4.*B1.*B2.*B5.^2-A2.*A3.*A5.^2.*B1.*B2.*B3.^2-A1.*A2.*A4.^2.*B2.*B3.^2.*B5+A1.*A2.*A5.^2.*B2.^2.*B3.*B4.*3.0-A1.*A3.*A4.^2.*B1.*B3.^2.*B5.*2.0+A1.*A3.*A5.^2.*B1.*B2.^2.*B5.*2.0+A1.*A4.*A5.^2.*B1.*B2.^2.*B4+A1.*A2.^2.*A3.*B2.*B4.*B5.^2+A1.*A2.^2.*A4.*B1.*B4.*B5.^2+A1.*A2.^2.*A4.*B2.*B3.*B5.^2.*2.0+A1.*A2.^2.*A5.*B1.*B3.*B5.^2.*2.0-A1.*A3.^2.*A5.*B1.*B3.*B4.^2.*2.0+A2.*A3.*A5.^2.*B1.*B2.^2.*B4.*2.0+A2.*A4.*A5.^2.*B1.*B2.^2.*B3-A2.*A3.^2.*A5.*B1.*B2.*B4.^2+A2.^2.*A3.*A4.*B1.*B2.*B5.^2.*3.0+A1.*A2.*A4.^2.*B2.^2.*B4.*B5.*2.0+A1.*A2.*A5.^2.*B1.^2.*B4.*B5.*5.0+A1.*A3.*A4.^2.*B2.^2.*B3.*B5+A1.*A3.*A5.^2.*B1.^2.*B3.*B5.*2.0-A1.*A4.*A5.^2.*B1.^2.*B2.*B5-A1.*A4.*A5.^2.*B1.^2.*B3.*B4.*5.0-A1.*A2.^2.*A4.*B2.*B4.^2.*B5-A1.*A2.^2.*A5.*B1.*B4.^2.*B5+A1.*A3.^2.*A5.*B1.*B3.^2.*B5.*4.0-A1.*A4.^2.*A5.*B1.*B2.^2.*B5-A2.*A3.*A5.^2.*B1.^2.*B2.*B5.*5.0+A2.*A3.*A5.^2.*B1.^2.*B3.*B4-A2.*A4.*A5.^2.*B1.^2.*B2.*B4-A2.*A4.^2.*A5.*B1.*B2.^2.*B4-A3.*A4.*A5.^2.*B1.^2.*B2.*B3-A1.^2.*A2.*A3.*B3.*B4.*B5.^2-A1.^2.*A2.*A4.*B2.*B4.*B5.^2-A1.^2.*A2.*A5.*B1.*B4.*B5.^2-A1.^2.*A2.*A5.*B2.*B3.*B5.^2.*5.0-A1.^2.*A3.*A4.*B1.*B4.*B5.^2.*5.0+A1.^2.*A3.*A4.*B2.*B3.*B5.^2+A1.^2.*A3.*A5.*B1.*B3.*B5.^2.*2.0+A1.^2.*A4.*A5.*B1.*B2.*B5.^2.*5.0+A2.^2.*A3.*A5.*B1.*B3.*B4.^2+A2.^2.*A4.*A5.*B1.*B2.*B4.^2.*2.0-A1.*A3.^2.*A4.*B2.^2.*B4.*B5-A1.*A3.^2.*A5.*B2.^2.*B3.*B5.*2.0+A1.*A4.^2.*A5.*B1.^2.*B3.*B5.*2.0+A2.*A3.*A4.^2.*B1.^2.*B4.*B5.*3.0+A2.*A4.^2.*A5.*B1.^2.*B2.*B5+A2.*A4.^2.*A5.*B1.^2.*B3.*B4.*2.0+A3.*A4.^2.*A5.*B1.^2.*B2.*B4+A1.^2.*A2.*A4.*B3.*B4.^2.*B5+A1.^2.*A2.*A5.*B2.*B4.^2.*B5+A1.^2.*A3.*A4.*B2.*B4.^2.*B5.*2.0+A1.^2.*A3.*A5.*B1.*B4.^2.*B5.*2.0+A1.^2.*A4.*A5.*B2.*B3.*B4.^2.*3.0-A2.^2.*A3.*A5.*B1.*B3.^2.*B5.*2.0-A2.^2.*A4.*A5.*B1.*B3.^2.*B4+A2.*A3.^2.*A5.*B1.^2.*B4.*B5+A1.^2.*A2.*A5.*B3.^2.*B4.*B5.*3.0-A1.^2.*A3.*A4.*B3.^2.*B4.*B5+A1.^2.*A4.*A5.*B2.*B3.^2.*B5+A3.^2.*A4.*A5.*B1.^2.*B2.*B5.*3.0-A3.^2.*A4.*A5.*B1.^2.*B3.*B4-A1.^2.*A4.*A5.*B2.^2.*B4.*B5.*5.0-A2.^2.*A4.*A5.*B1.^2.*B4.*B5.*5.0-A2.^2.*A5.^2.*B1.*B2.*B3.*B4.*3.0-A1.^2.*A5.^2.*B1.*B2.*B4.*B5.*8.0-A2.^2.*A4.^2.*B1.*B2.*B4.*B5.*2.0-A1.^2.*A4.^2.*B2.*B3.*B4.*B5.*3.0-A1.*A2.*A3.*A4.*B2.^2.*x2.^2.*3.0-A1.*A2.*A4.*A5.*B1.^2.*x2.^2.*8.0-A1.*A2.*A4.^2.*B1.*B2.*x2.^2.*5.0+A1.*A2.*A3.^2.*B1.*B4.*x2.^2.*3.0-A1.*A2.*A3.^2.*B2.*B3.*x2.^2+A1.*A3.^2.*A4.*B1.*B2.*x2.^2+A1.*A2.*A4.^2.*B2.*B3.^2.*x2+A1.*A3.*A4.^2.*B1.*B3.^2.*x2.*2.0-A1.*A3.*A5.^2.*B1.*B2.^2.*x2.*2.0-A1.*A2.^2.*A3.*B1.*B5.*x2.^2.*1.2e1+A1.*A2.^2.*A3.*B1.*B5.^2.*x2.*1.2e1+A1.*A2.^2.*A3.*B2.*B4.*x2.^2+A1.*A2.^2.*A4.*B1.*B4.*x2.^2+A1.*A2.^2.*A4.*B2.*B3.*x2.^2.*2.0+A1.*A2.^2.*A5.*B1.*B3.*x2.^2.*2.0+A2.^2.*A3.*A4.*B1.*B2.*x2.^2.*3.0-A1.*A2.*A4.^2.*B2.^2.*B4.*x2.*2.0-A1.*A2.*A5.^2.*B1.^2.*B4.*x2.*5.0-A1.*A3.*A4.^2.*B2.^2.*B3.*x2-A1.*A3.*A5.^2.*B1.^2.*B3.*x2.*2.0+A1.*A4.*A5.^2.*B1.^2.*B2.*x2+A1.*A2.^2.*A4.*B2.*B4.^2.*x2+A1.*A2.^2.*A5.*B1.*B4.^2.*x2-A1.*A3.^2.*A5.*B1.*B3.^2.*x2.*4.0+A1.*A4.^2.*A5.*B1.*B2.^2.*x2+A2.*A3.*A5.^2.*B1.^2.*B2.*x2.*5.0+A1.^2.*A2.*A3.*B2.*B5.*x2.^2.*9.0-A1.^2.*A2.*A3.*B2.*B5.^2.*x2.*9.0-A1.^2.*A2.*A3.*B3.*B4.*x2.^2+A1.^2.*A2.*A4.*B1.*B5.*x2.^2.*1.2e1-A1.^2.*A2.*A4.*B1.*B5.^2.*x2.*1.2e1-A1.^2.*A2.*A4.*B2.*B4.*x2.^2-A1.^2.*A2.*A5.*B1.*B4.*x2.^2-A1.^2.*A2.*A5.*B2.*B3.*x2.^2.*5.0-A1.^2.*A3.*A4.*B1.*B4.*x2.^2.*5.0+A1.^2.*A3.*A4.*B2.*B3.*x2.^2+A1.^2.*A3.*A5.*B1.*B3.*x2.^2.*2.0+A1.^2.*A4.*A5.*B1.*B2.*x2.^2.*5.0-A1.*A3.*A4.^2.*B1.^2.*B5.*x2.*8.0+A1.*A3.^2.*A4.*B2.^2.*B4.*x2+A1.*A3.^2.*A5.*B2.^2.*B3.*x2.*2.0-A1.*A4.^2.*A5.*B1.^2.*B3.*x2.*2.0-A2.*A3.*A4.^2.*B1.^2.*B4.*x2.*3.0-A2.*A4.^2.*A5.*B1.^2.*B2.*x2-A1.^2.*A2.*A4.*B3.*B4.^2.*x2-A1.^2.*A2.*A5.*B2.*B4.^2.*x2-A1.^2.*A3.*A4.*B2.*B4.^2.*x2.*2.0-A1.^2.*A3.*A5.*B1.*B4.^2.*x2.*2.0+A2.^2.*A3.*A5.*B1.*B3.^2.*x2.*2.0-A1.*A2.^2.*A5.*B2.^2.*B5.*x2.*6.0+A1.*A3.^2.*A5.*B1.^2.*B5.*x2.*8.0+A2.*A3.^2.*A4.*B1.^2.*B5.*x2.*8.0-A2.*A3.^2.*A5.*B1.^2.*B4.*x2+A1.^2.*A2.*A4.*B3.^2.*B5.*x2.*4.0-A1.^2.*A2.*A5.*B3.^2.*B4.*x2.*3.0+A1.^2.*A3.*A4.*B3.^2.*B4.*x2-A1.^2.*A4.*A5.*B2.*B3.^2.*x2-A3.^2.*A4.*A5.*B1.^2.*B2.*x2.*3.0+A1.^2.*A3.*A5.*B2.^2.*B5.*x2.*6.0+A1.^2.*A4.*A5.*B2.^2.*B4.*x2.*5.0-A2.^2.*A3.*A5.*B1.^2.*B5.*x2.*8.0+A2.^2.*A4.*A5.*B1.^2.*B4.*x2.*5.0+A1.^2.*A5.^2.*B1.*B2.*B4.*x2.*8.0+A2.^2.*A4.^2.*B1.*B2.*B4.*x2.*2.0+A1.^2.*A4.^2.*B1.*B3.*B5.*x2.*6.0+A1.^2.*A4.^2.*B2.*B3.*B4.*x2.*3.0-A2.^2.*A3.^2.*B1.*B3.*B5.*x2.*2.0+A1.^2.*A3.^2.*B2.*B4.*B5.*x2.*4.0+A1.*A2.*A3.*A4.*B1.*B3.*B5.^2.*4.0+A1.*A2.*A3.*A5.*B1.*B2.*B5.^2.*2.0-A1.*A2.*A3.*A4.*B1.*B4.^2.*B5.*3.0-A1.*A2.*A3.*A5.*B2.*B3.*B4.^2-A1.*A2.*A4.*A5.*B1.*B3.*B4.^2+A1.*A3.*A4.*A5.*B1.*B2.*B4.^2+A1.*A2.*A3.*A5.*B2.*B3.^2.*B5.*2.0+A1.*A2.*A4.*A5.*B2.*B3.^2.*B4+A1.*A3.*A4.*A5.*B1.*B3.^2.*B4.*2.0+A1.*A2.*A3.*A5.*B2.^2.*B4.*B5-A1.*A2.*A4.*A5.*B2.^2.*B3.*B5-A1.*A3.*A4.*A5.*B2.^2.*B3.*B4-A2.*A3.*A4.*A5.*B1.*B2.^2.*B5.*3.0+A1.*A3.*A4.*A5.*B1.^2.*B4.*B5.*2.0+A2.*A3.*A4.*A5.*B1.^2.*B3.*B5.*4.0+A1.*A2.*A5.^2.*B1.*B2.*B3.*B5.*2.0+A1.*A3.*A5.^2.*B1.*B2.*B3.*B4.*4.0+A1.*A2.*A4.^2.*B1.*B3.*B4.*B5-A1.*A3.*A4.^2.*B1.*B2.*B4.*B5-A1.*A4.^2.*A5.*B1.*B2.*B3.*B4.*3.0-A2.*A3.*A4.^2.*B1.*B2.*B3.*B5+A1.*A3.^2.*A4.*B1.*B3.*B4.*B5.*2.0+A2.*A3.^2.*A4.*B1.*B2.*B4.*B5+A2.*A3.^2.*A5.*B1.*B2.*B3.*B5.*2.0-A1.*A2.^2.*A5.*B2.*B3.*B4.*B5.*3.0-A2.^2.*A3.*A4.*B1.*B3.*B4.*B5-A2.^2.*A3.*A5.*B1.*B2.*B4.*B5+A2.^2.*A4.*A5.*B1.*B2.*B3.*B5+A1.^2.*A3.*A5.*B2.*B3.*B4.*B5.*4.0+A1.^2.*A4.*A5.*B1.*B3.*B4.*B5.*2.0+A1.*A2.*A3.*A4.*B1.*B3.*x2.^2.*4.0+A1.*A2.*A3.*A5.*B1.*B2.*x2.^2.*2.0+A1.*A2.*A3.*A4.*B1.*B4.^2.*x2.*3.0-A1.*A2.*A3.*A5.*B2.*B3.^2.*x2.*2.0+A1.*A2.*A3.*A4.*B2.^2.*B5.*x2.*6.0-A1.*A2.*A3.*A5.*B2.^2.*B4.*x2+A1.*A2.*A4.*A5.*B2.^2.*B3.*x2+A2.*A3.*A4.*A5.*B1.*B2.^2.*x2.*3.0+A1.*A2.*A4.*A5.*B1.^2.*B5.*x2.*1.6e1-A1.*A3.*A4.*A5.*B1.^2.*B4.*x2.*2.0-A2.*A3.*A4.*A5.*B1.^2.*B3.*x2.*4.0-A1.*A2.*A5.^2.*B1.*B2.*B3.*x2.*2.0+A1.*A2.*A4.^2.*B1.*B2.*B5.*x2.*1.0e1-A1.*A2.*A4.^2.*B1.*B3.*B4.*x2+A1.*A3.*A4.^2.*B1.*B2.*B4.*x2+A2.*A3.*A4.^2.*B1.*B2.*B3.*x2-A1.*A2.*A3.^2.*B1.*B4.*B5.*x2.*6.0+A1.*A2.*A3.^2.*B2.*B3.*B5.*x2.*2.0-A1.*A3.^2.*A4.*B1.*B2.*B5.*x2.*2.0-A1.*A3.^2.*A4.*B1.*B3.*B4.*x2.*2.0-A2.*A3.^2.*A4.*B1.*B2.*B4.*x2-A2.*A3.^2.*A5.*B1.*B2.*B3.*x2.*2.0-A1.*A2.^2.*A3.*B2.*B4.*B5.*x2.*2.0-A1.*A2.^2.*A4.*B1.*B4.*B5.*x2.*2.0-A1.*A2.^2.*A4.*B2.*B3.*B5.*x2.*4.0-A1.*A2.^2.*A5.*B1.*B3.*B5.*x2.*4.0+A1.*A2.^2.*A5.*B2.*B3.*B4.*x2.*3.0-A2.^2.*A3.*A4.*B1.*B2.*B5.*x2.*6.0+A2.^2.*A3.*A4.*B1.*B3.*B4.*x2+A2.^2.*A3.*A5.*B1.*B2.*B4.*x2-A2.^2.*A4.*A5.*B1.*B2.*B3.*x2+A1.^2.*A2.*A3.*B3.*B4.*B5.*x2.*2.0+A1.^2.*A2.*A4.*B2.*B4.*B5.*x2.*2.0+A1.^2.*A2.*A5.*B1.*B4.*B5.*x2.*2.0+A1.^2.*A2.*A5.*B2.*B3.*B5.*x2.*1.0e1+A1.^2.*A3.*A4.*B1.*B4.*B5.*x2.*1.0e1-A1.^2.*A3.*A4.*B2.*B3.*B5.*x2.*2.0-A1.^2.*A3.*A5.*B1.*B3.*B5.*x2.*4.0-A1.^2.*A3.*A5.*B2.*B3.*B4.*x2.*4.0-A1.^2.*A4.*A5.*B1.*B2.*B5.*x2.*1.0e1-A1.^2.*A4.*A5.*B1.*B3.*B4.*x2.*2.0+A1.*A2.*A3.*A4.*B2.*B3.*B4.*B5-A1.*A2.*A3.*A5.*B1.*B3.*B4.*B5.*8.0+A1.*A2.*A4.*A5.*B1.*B2.*B4.*B5.*1.0e1-A1.*A3.*A4.*A5.*B1.*B2.*B3.*B5.*8.0+A2.*A3.*A4.*A5.*B1.*B2.*B3.*B4-A1.*A2.*A3.*A4.*B1.*B3.*B5.*x2.*8.0-A1.*A2.*A3.*A4.*B2.*B3.*B4.*x2-A1.*A2.*A3.*A5.*B1.*B2.*B5.*x2.*4.0+A1.*A2.*A3.*A5.*B1.*B3.*B4.*x2.*8.0-A1.*A2.*A4.*A5.*B1.*B2.*B4.*x2.*1.0e1+A1.*A3.*A4.*A5.*B1.*B2.*B3.*x2.*8.0;
