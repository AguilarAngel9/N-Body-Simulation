%% Entregable 4. Problema de 3 cuerpos.
function ThreeBodyProblem()
    %% Inicialización
    clear
    clc
    %Constantes
    G = 6.67*10^-11;
    %Unidad astronomica
    UA = 1.496*10^11;
    
    %Masas
    M_Sol = 1.989*10^30;
    M_Jup = 1.898*10^27;
    M_Ast = 6*10^12;
    
    %Caso a realizar
    caso = "Atrapamiento";
    
    switch caso
        %Se declaran las condiciones iniciales dependiendo del caso. 
        case 'Orbita'
            ro1 = [0,0,0]; % Sol
            ro2 = [740.68*10^9,0,0]; % Jup
            ro3 = [740.68*10^9,-2.2*10^9,0]; % Ast
            
            vo1 = [0,0,0]; % Sol 
            vo2 = [0,13.72*10^3,0]; % Jup
            vo3 = [-10000,13.72*10^3,0]; % Ast
            
            %Tiempo y iteraciones
            h = 19900; %Paso de tiempo
            t = 20500; %Iteraciones
            
        case 'Atrapamiento'

            ro1 = [0,0,0]; % Sol
            ro2 = [740.68*10^9,0,0]; % Jup
            ro3 = [6.0325*UA, 5.63*UA, 0]; % Ast
            
            vo1 = [0,0,0]; % Sol
            vo2 = [0,13.72*10^3,0]; % Jup
            vo3 = [-7000,0,0]; % Ast
            
            %Tiempo y iteraciones
            h = 16000; %Paso de tiempo
            t = 25500; %Iteraciones
            
        case 'Orbita en z'
            ro1 = [0,0,0]; % Sol
            ro2 = [740.68*10^9,0,0]; % Jup
            ro3 = [740.68*10^9,-2.2*10^9,0]; % Ast
            
            vo1 = [0,0,0]; % Sol
            vo2 = [0,13.72*10^3,0]; % Jup
            vo3 = [500,13.82*10^3,-500]; % Ast
            
            %Tiempo y iteraciones
            h = 2500;
            t = 150000;
    end

    %Animación
    key = 0;
    
    %Centros de masas
    MT = M_Sol + M_Jup + M_Ast;
    RCM = (ro1*M_Sol + ro2*M_Jup + ro3*M_Ast)/MT;
    
    %% Vectores a utilizar
    
    %Regresan los vectores distancia escalar, vectorial dependiendo de 
    % los objetos y la posición respescto al origen
    
    %Vectores posición del Sol
    [r_Sol, DSol_Jup, DSol_Ast, rJup_Sol, rAst_Sol] = Distancias(ro1,ro2,ro3,t);
    
    %Vectores posición de Jupiter
    [r_Jup, ~, DJup_Ast, rSol_Jup, rAst_Jup] = Distancias(ro2,ro1,ro3,t);
    
    %Vectores posición del asteroide
    [r_Ast, ~, ~, rSol_Ast, rJup_Ast] = Distancias(ro3,ro1,ro2,t);
    
    
    %Regresan vectores de velocidad de cada objeto.
    
    %Vector velocidad del sol
    [v_Sol, v_half_Sol] = velocidades(vo1,t);
    
    %Vector velocidad de jupiter
    [v_Jup, v_half_Jup] = velocidades(vo2,t);
    
    %Vector velocidad del asteroide
    [v_Ast, v_half_Ast] = velocidades(vo3,t);
    
    %Se declaran las fuerzas y su primer valor en t=0.
    %Fuerzas y su calculo inicial
    F_Sol = zeros(t,3);
    F_Jup = zeros(t,3);
    F_Ast = zeros(t,3);
    
    [F_Sol] = fuerzas(F_Sol, G, M_Sol, M_Jup, M_Ast, DSol_Jup(1,:), DSol_Ast(1,:), rJup_Sol, rAst_Sol, 1);
    [F_Jup] = fuerzas(F_Jup, G, M_Jup, M_Sol, M_Ast, DSol_Jup(1,:), DJup_Ast(1,:), rSol_Jup, rAst_Jup, 1);
    [F_Ast] = fuerzas(F_Ast, G, M_Ast, M_Sol, M_Jup, DSol_Ast(1,:), DJup_Ast(1,:), rSol_Ast, rJup_Ast, 1);
    
    
    %% Velocity Verlet
    
    %Inicializan contadores para la animación y frenado si llega a ser el
    %caso de atrapamiento.
    
    %Animación
    j = 1;
    if key ==1
        %Despliega una figura para la animación si se solicita con key.
        figure(2)
    end
    
    %Frenado
    cont = 1;
    
    %Inicio del bucle para calculos de verlet.
    for i=2:t
        %Velocidades medias
        v_half_Sol(i-1,:) = v_Sol(i-1,:) + 0.5* h *F_Sol(i-1,:)/M_Sol;
        v_half_Jup(i-1,:) = v_Jup(i-1,:) + 0.5* h *F_Jup(i-1,:)/M_Jup;
        v_half_Ast(i-1,:) = v_Ast(i-1,:) + 0.5* h *F_Ast(i-1,:)/M_Ast;
        
        %Posiciones nuevas
        r_Sol(i,:) = r_Sol(i-1,:) + h * v_half_Sol(i-1,:);
        r_Jup(i,:) = r_Jup(i-1,:) + h * v_half_Jup(i-1,:);
        r_Ast(i,:) = r_Ast(i-1,:) + h * v_half_Ast(i-1,:);
        
        %Nueva posición del centro de masas
        RCM(i,:) = (r_Sol(i,:)*M_Sol + r_Jup(i,:)*M_Jup + r_Ast(i,:)*M_Ast)/MT;
        
        %Vectores unitarios y distancia escalar y vectorial para cada 
        %cuerpo respecto a los demas cuerpos.
        %Se llama a la función unitarios que regresa todos esos calculos.
        
        %Sol
        [rJup_Sol, rAst_Sol, DSol_Jup, DSol_Ast] = unitarios(r_Sol(i,:), r_Jup(i,:), r_Ast(i,:), DSol_Jup, DSol_Ast, i);
        %Jup
        [rSol_Jup, rAst_Jup, ~, DJup_Ast] = unitarios(r_Jup(i,:), r_Sol(i,:), r_Ast(i,:), DSol_Jup, DJup_Ast, i);
        
        %Ast
        %Como los vectores posición respecto a los cuerpos son iguales pero
        %de signo contrario, solo se multiplica por -1 los del asteroide.
        rSol_Ast = -rAst_Sol;
        rJup_Ast = -rAst_Jup;
        
        %Aquí se ejecuta el frenado del asteroide si se cumple que la
        %distancia entre jupiter y el asteroide son menores a una distancia
        %este se ejecutara.
        if ((caso == "Atrapamiento") && (DJup_Ast(i) < 1*10^9) && (cont == 1))
            %El contador se hace cero ya que solo frena una sola vez.
            cont = 0;
            %Se calcula el signo del vector velocidad 
            signo = v_Ast(i-1,[1,2])./abs(v_Ast(i-1,[1,2]));
            
            %Se realiza el frenado, el signo es para obtener una suma en la
            %dirección opuesta y solo es un quinto de su velocidad la que
            %se resta.
            v_Ast(i,[1,2]) = v_Ast(i-1,[1,2]) + signo.*v_Ast(i-1,[1,2])/5;

            %Se calcula la fuerza del asteroide
            F_Ast = fuerzas(F_Ast, G, M_Ast, M_Sol, M_Jup, DSol_Ast(i), DJup_Ast(i), rSol_Ast, rJup_Ast, i);
            
        else
            %Si no se cumplen las condiciones, ya sea el caso o la
            %distancia, entonces se calcula la Fuerza y velocidad del
            %asteroide de una manera normal, sin el frenado.
            F_Ast = fuerzas(F_Ast, G, M_Ast, M_Sol, M_Jup, DSol_Ast(i), DJup_Ast(i), rSol_Ast, rJup_Ast, i);
            v_Ast(i,:) = v_half_Ast(i-1,:) + 0.5*h*F_Ast(i,:)/M_Ast;
        end
        
        %Fuerzas
        F_Sol = fuerzas(F_Sol, G, M_Sol, M_Jup, M_Ast, DSol_Jup(i), DSol_Ast(i), rJup_Sol, rAst_Sol, i);
        F_Jup = fuerzas(F_Jup, G, M_Jup, M_Sol, M_Ast, DSol_Jup(i), DJup_Ast(i), rSol_Jup, rAst_Jup, i);
        
        %Velocidades
        v_Sol(i,:) = v_half_Sol(i-1,:) + 0.5*h*F_Sol(i,:)/M_Sol;
        v_Jup(i,:) = v_half_Jup(i-1,:) + 0.5*h*F_Jup(i,:)/M_Jup;
        
        %Si se requiere observar la animación, aquí se ejecuta y grafica.
        if j == 300 && key == 1
            animacion(r_Sol, r_Jup, r_Ast);
            j = 0;
        end
        j = j + 1;
    end
    %% Graficas resultantes.
    figure(1)
    %Se grafica el sol
    scatter3(r_Sol(:,1), r_Sol(:,2), r_Sol(:,3), 'filled', 'r', 'DisplayName', 'Sol')
    hold on
    %Se grafican las trayectorias de jupiter y el asteroide.
    plot3(r_Jup(:,1), r_Jup(:,2), r_Jup(:,3), 'color', '#ad4234', 'DisplayName', 'Jupiter')
    plot3(r_Ast(:,1), r_Ast(:,2), r_Ast(:,3), 'color', '#6c45a3', 'DisplayName', 'Asteroide')
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend
    grid off
    axis equal
    
    %Se imprime la distancia más corta entre el asteroide y jupiter que
    %existio.
    fprintf('Minima, %d \n',min(DJup_Ast))
end


function [rp, s1_2, s1_3, runit_1, runit_2] = Distancias(ro1,ro2,ro3,t)
    %Esta función declara los vectores posición, distancia escalar, y los
    %unitarios del objeto que se solicite respecto a los otros 2.
    rp = zeros(t,3);
    s1_2 = zeros(t,1);
    s1_3 = zeros(t,1);
    
    rp(1,:) = ro1;
    
    %Cabe destacar que cuerpo 1, 2, 3, van cambiando dependendiendo del
    %tipo de entrada que se le da a la función
    %El 1 puede ser sol, o jupiter o el asteroide.
    
    %Vector distancia del cuerpo 2 al 1 y del cuerpo 3 al 1.
    r1 = ro1-ro2;
    r2 = ro1-ro3;
    
    %Distancia escalar del cuerpo 1 al 2 y al 3.
    s1_2(1) = norm(r1);
    s1_3(1) = norm(r2);
    
    %Unitarios del vector distancia.
    runit_1 = r1/s1_2(1);
    runit_2 = r2/s1_3(1);
end

function [v, v_half] = velocidades(vo,t)
    %Función que declara los vectores velocidad y mitad de velocidad con la
    %condición inicial.
    v = zeros(t,3);
    v_half = zeros(t,3);
    v(1,:) = vo;
end

function [F] = fuerzas(F, G, M1, M2, M3, s1, s2, runit1, runit2, i)
    %Hace los calculos de las fuerzas, entre el cuerpo 1 con el 2 y el 3
    %Cabe destacar que el cuerpo 1,2,3 pueden rotar, dependen de la entrada
    %que le des a la función
    %Cuerpo 1 puede ser sol, jup, asteroide.
    
    F12 = -((G*M1*M2)/(s1^2))*runit1;
    F13 = -((G*M1*M3)/(s2^2))*runit2;
    
    %Se realiza la suma de fuerzas y es lo que regresa la función.
    F(i,:) = F12 + F13;
end

function [runit_1, runit_2, s1_2, s1_3] = unitarios(r1, r2, r3, s1_2, s1_3, i)
    %Hace el calculo de los nuevos vectores posición y unitarios para cada
    %intervalo i.
    rp1 = r1-r2;
    rp2 = r1-r3;
    
    s1_2(i) = norm(rp1);
    s1_3(i) = norm(rp2);
    
    runit_1 = rp1/s1_2(i);
    runit_2 = rp2/s1_3(i);
    
end

function animacion(r_sol, r_jup, r_ast)
    %Función que realiza la animación de la evolución de las trayectorias,
    %para fines visuales de las
    
    scatter3(r_sol(:,1), r_sol(:,2), r_sol(:,3), 2, 'filled')
    hold on
    scatter3(r_jup(:,1), r_jup(:,2), r_jup(:,3), 1, 'filled', 'r')
    scatter3(r_ast(:,1), r_ast(:,2), r_ast(:,3), 1, 'filled', 'b')
    pause(0.01)
end
