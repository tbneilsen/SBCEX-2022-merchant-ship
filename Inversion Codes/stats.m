% /****************************************************************
%  *                                                              *
%  *     COPYRIGHT © 2016                                         *
%  *     Knobles Scientific and Analysis, LLC                     *                                      *
%  *     All rights reserved                                      *
%  *                                                              *
%  * Redistribution of source or binary forms is not permitted    *
%  * without the permission of KSA, LLC                           *
%  * THIS SOFTWARE IS PROVIDED BY KSA,LLC “AS IS” AND ANY EXPRESS *
%  * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    *
%  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A      *
%  * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL KSA,LLC *
%  * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
%  * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT      *
%  * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
%  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)     *
%  * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN    *
%  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR *
%  * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS         *
%  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
%  *                                                              *
%  ****************************************************************/


function [Stat] = stats(WI,PWI,nvar)

%Loop throu each parameter and obtain the mean and standard deviation
%Note that both must be weighted by the probability distribution
for i = 1:nvar

    %The sum(P(PWI)) should equal 1.  Obtain the 
    %normalization constant, N, after the fact.
    N(i) = sum(PWI(i,:));
    %Normalize by dividing by N.
    P(i,:) = PWI(i,:)./N(i);
    %Compute the mean
    Stat(i,1) = ( WI(i,:)*P(i,:)' );
    %Compute the standard deviation
    Stat(i,2) = sqrt( (WI(i,:) - Stat(i,1) ).^2*P(i,:)' );
    
end

%Transform stats into a full matrix instead of a sparse matrix
Stat = full(Stat);